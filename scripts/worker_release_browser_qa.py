#!/usr/bin/env python3
"""Prove a complete content-addressed worker release with the existing HTTP cache.

Old fixed URLs are armed with invalid, long-lived-cache sentinels. They are fetched
inertly before fresh navigation to verify that the guard is active, then any use
by the application fails the run. No HTTP cache clear, cache-disable or profile
replacement is used. All fresh startup Runtime/console/Worker errors are retained.
"""
import base64
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urljoin, urlparse

from browser_qa import QA_DIST, QA_ORIGIN, Browser
import desktop_mvp_browser_qa as mvp

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = Path('/tmp/collagene2-worker-release-qa')
ORIGIN = QA_ORIGIN
NAMES = ('analysis-worker_loader.js', 'analysis-worker.js', 'analysis-worker_bg.wasm')
OLD_URLS = {ORIGIN + '/' + name for name in NAMES}
PATTERNS = [{'urlPattern': url + '*', 'requestStage': 'Request'} for url in sorted(OLD_URLS)]
SENTINEL = 'STALE_FIXED_WORKER_SENTINEL'

# Observe transport identity only. This does not decode any scientific payload.
PROBE = r'''
(() => {
    const Original = window.Worker;
    window.__releaseOriginalWorker = Original;
    window.__releaseWorkers = [];
    const summary = data => {
        const bytes = data instanceof Uint8Array ? data : new Uint8Array(data);
        if (bytes.length < 42) throw Error('Short worker envelope');
        let at = 42;
        const integer = () => {
            let value = 0, scale = 1;
            for (let i = 0; i < 8; ++i) {
                if (at >= bytes.length) throw Error('Truncated envelope identity');
                const next = bytes[at++]; value += (next & 127) * scale;
                if ((next & 128) === 0) return value;
                scale *= 128;
            }
            throw Error('Oversized envelope identity');
        };
        const id = integer(), revisionLength = integer();
        if (revisionLength > 4096 || at + revisionLength >= bytes.length) throw Error('Invalid revision identity');
        const revision = new TextDecoder().decode(bytes.slice(at, at + revisionLength));
        at += revisionLength;
        return {magic: String.fromCharCode(...bytes.slice(0, 4)), version: bytes[4] + bytes[5] * 256,
            id, revision, tag: integer(), bytes: bytes.length};
    };
    window.Worker = new Proxy(Original, {
        construct(Target, args) {
            const worker = Reflect.construct(Target, args);
            const entry = {url: new URL(String(args[0]), location.href).href, ready: false,
                outgoing: [], incoming: [], controls: [], errors: []};
            window.__releaseWorkers.push(entry);
            const post = worker.postMessage;
            worker.postMessage = function(data, ...rest) {
                entry.outgoing.push(summary(data));
                return post.call(this, data, ...rest);
            };
            worker.addEventListener('message', event => {
                if (typeof event.data === 'string') {
                    if (event.data === 'ready') entry.ready = true;
                    else entry.controls.push(event.data);
                } else entry.incoming.push(summary(event.data));
            });
            worker.addEventListener('error', error => entry.errors.push(error.message || 'Worker error'));
            return worker;
        }
    });
})();
'''


class TransportTape:
    """Pass every CDP message through; retain child-worker requests separately.

    Browser.call still receives all original Runtime/Log events. Flattened worker
    sessions enable Network before resuming new workers. Chrome104 has no child
    Fetch domain; child paths are checked from actual Network events, while the
    page Fetch guard seeds and rejects legacy resources.
    """
    def __init__(self, socket):
        self.socket = socket
        self.serial = 1_000_000_000
        self.commands = {}
        self.sessions = {}
        self.events = []
        self.legacy_hits = []
        self.control_errors = []

    def __getattr__(self, name):
        return getattr(self.socket, name)

    def command(self, method, params, session=None):
        self.serial += 1
        packet = {'id': self.serial, 'method': method, 'params': params}
        if session:
            packet['sessionId'] = session
        self.commands[self.serial] = (method, session)
        self.socket.send(json.dumps(packet))

    def recv(self, *args, **kwargs):
        raw = self.socket.recv(*args, **kwargs)
        event = json.loads(raw)
        method, params = event.get('method'), event.get('params', {})
        session = event.get('sessionId')
        if method == 'Target.attachedToTarget':
            attached = params['sessionId']
            self.sessions[attached] = params.get('targetInfo', {})
            self.command('Network.enable', {}, attached)
            self.command('Runtime.runIfWaitingForDebugger', {}, attached)
        if method == 'Fetch.requestPaused':
            request = params['request']
            url = request['url']
            if url.split('?', 1)[0] in OLD_URLS:
                self.legacy_hits.append({'url': url, 'session': session})
                body = SENTINEL if urlparse(url).path.endswith('.wasm') else 'throw new Error(' + json.dumps(SENTINEL) + ');'
                self.command('Fetch.fulfillRequest', {
                    'requestId': params['requestId'], 'responseCode': 200,
                    'responseHeaders': [
                        {'name': 'Content-Type', 'value': 'application/wasm' if urlparse(url).path.endswith('.wasm') else 'application/javascript'},
                        {'name': 'Cache-Control', 'value': 'public, max-age=31536000'},
                    ],
                    'body': base64.b64encode(body.encode()).decode(),
                }, session)
            else:
                self.command('Fetch.continueRequest', {'requestId': params['requestId']}, session)
        if method in ('Network.requestWillBeSent', 'Network.responseReceived', 'Network.loadingFailed'):
            payload = params.get('request', params.get('response', {}))
            self.events.append({'method': method, 'url': payload.get('url'), 'session': session,
                'target': self.sessions.get(session, {}).get('url'),
                'from_disk_cache': payload.get('fromDiskCache', False),
                'error': params.get('errorText')})
            if len(self.events) > 20_000:
                raise AssertionError('Unexpected unbounded QA network-event stream')
        if event.get('id') in self.commands:
            command = self.commands.pop(event['id'])
            if 'error' in event:
                self.control_errors.append({'command': command, 'error': event['error']})
        return raw


def assert_worker_release(browser, tape, name, isolated_tags=(19,)):
    meta = browser.evaluate('document.querySelector("meta[name=collagene-worker]")?.content')
    assert meta, 'Worker release meta is missing'
    absolute = urljoin(ORIGIN + '/', meta)
    assert urlparse(absolute).netloc == urlparse(ORIGIN).netloc
    assert re.fullmatch(r'.*/workers/[0-9a-f]{64}/analysis-worker_loader\.js', absolute), absolute
    directory = absolute.rsplit('/', 1)[0] + '/'
    expected = {directory + filename for filename in NAMES}
    entries = browser.evaluate('window.__releaseWorkers')
    assert entries, 'No real workers were constructed'
    assert all(entry['url'] == absolute for entry in entries), entries
    assert all(not entry['controls'] and not entry['errors'] for entry in entries), entries
    # Source changes can terminate an unused worker before its ready handshake.
    # Every worker that actually exchanged an envelope must have become ready.
    assert all(entry['ready'] for entry in entries if entry['outgoing'] or entry['incoming']), entries
    for entry in entries:
        assert all(item['magic'] == 'CGQ\x00' and item['version'] == 22 for item in entry['outgoing']), entry
        assert all(item['magic'] == 'CGR\x00' and item['version'] == 22 for item in entry['incoming']), entry
        for reply in entry['incoming']:
            assert any(request['id'] == reply['id'] and request['revision'] == reply['revision'] for request in entry['outgoing']), entry
    # Stable appended wire identities:19 SaveWorkspace /20 LoadWorkspace (isolated),32 PrepareTracks,
    #33 CachedViewTracks (same retained worker),70 InspectCatalogRestriction.
    def answered(entry, tag):
        return entry['ready'] and any(request['tag'] == tag and any(
            reply['id'] == request['id'] and reply['revision'] == request['revision']
            for reply in entry['incoming']) for request in entry['outgoing'])
    assert any(any(answered(entry, tag) for tag in isolated_tags) for entry in entries), entries
    assert any(answered(entry, 32) and answered(entry, 33) for entry in entries), entries
    assert any(answered(entry, 70) for entry in entries), entries
    worker_events = [event for event in tape.events if event['url'] and urlparse(event['url']).path.rsplit('/', 1)[-1] in NAMES]
    assert worker_events, 'No worker network evidence was retained'
    assert all(event['url'] in expected for event in worker_events), worker_events
    observed = {event['url'] for event in worker_events}
    assert directory + 'analysis-worker.js' in observed and directory + 'analysis-worker_bg.wasm' in observed, worker_events
    assert not tape.legacy_hits, tape.legacy_hits
    assert not tape.control_errors, tape.control_errors
    assert not browser.exceptions and not browser.console_errors, (browser.exceptions, browser.console_errors)
    (DIRECTORY / (name + '.json')).write_text(json.dumps({'meta': meta, 'workers': entries,
        'network': tape.events, 'legacy_requests': tape.legacy_hits, 'cdp_control_errors': tape.control_errors}, indent=2))
    return meta


def inspect_catalog(browser):
    mvp.click(browser, 'Digest')
    mvp.set_enzymes(browser, ['EcoRI'])
    mvp.click(browser, 'Run virtual digest')
    browser.wait('!!document.querySelector(".digest-results")', timeout=60)
    artifact = json.loads(mvp.download(browser, 'Download digest evidence', 'virtual-digest.json'))
    report = artifact['report']
    assert report['scan_complete'] and report['products_determined']
    assert [(site['enzyme'], site['recognition_start'], site['top_cut'], site['bottom_cut']) for site in report['sites']] == [('EcoRI', 4, 5, 9)]
    assert report['ref_top_lengths'] == [11, 5]
    assert artifact['catalog_engine_version'] == 'catalog-digest-1'
    return artifact


def main():
    DIRECTORY.mkdir(exist_ok=True)
    mvp.DIRECTORY = DIRECTORY
    browser = Browser()
    browser.socket = tape = TransportTape(browser.socket)
    browser.call('Browser.setDownloadBehavior', {'behavior': 'allow', 'downloadPath': str(DIRECTORY)})
    browser.call('Network.enable')
    browser.call('Fetch.enable', {'patterns': PATTERNS})
    browser.call('Target.setAutoAttach', {'autoAttach': True, 'waitForDebuggerOnStart': True, 'flatten': True})
    hook = browser.call('Page.addScriptToEvaluateOnNewDocument', {'source': PROBE})['identifier']
    try:
        # Inert fetches verify the invalid old-resource guard and seed old-path
        # cache entries where Chrome caches fulfilled requests. The guard remains
        # armed even when an existing cache entry bypasses ordinary revalidation.
        for url in sorted(OLD_URLS):
            text = browser.evaluate('(async()=>await(await fetch(' + json.dumps(url) + ',{cache:"reload"})).text())()')
            assert SENTINEL in text
        assert len(tape.legacy_hits) >= 3
        browser.call('Storage.clearDataForOrigin', {'origin': ORIGIN, 'storageTypes': 'indexeddb'})
        # Clear historical Runtime.enable replay BEFORE starting the fresh app.
        # Never clear fresh startup exceptions after navigation/readiness.
        browser.exceptions.clear()
        browser.console_errors.clear()
        tape.events.clear()
        tape.legacy_hits.clear()
        browser.call('Page.navigate', {'url': ORIGIN})
        browser.wait('document.body?.innerText.includes("Saved in this browser")')
        browser.wait('document.querySelector(".track-status")?.textContent.includes("Visible tracks")')
        sequence = 'AAAAGAATTCTTTCCC'
        path = DIRECTORY / 'WORKER_RELEASE.gb'
        path.write_text(mvp.record('WORKER_RELEASE', sequence, False))
        mvp.upload(browser, path)
        browser.wait('document.querySelector("h1")?.textContent==="WORKER_RELEASE"')
        browser.wait('document.querySelector(".track-status")?.textContent.includes("Visible tracks")')
        inspect_catalog(browser)
        browser.wait('document.body?.innerText.includes("Saved in this browser")')
        meta = assert_worker_release(browser, tape, 'first-startup')
        # Browser-fetched staged assets must equal this exact distribution's bytes.
        asset_directory = urljoin(ORIGIN + '/', meta).rsplit('/', 1)[0] + '/'
        for filename in NAMES:
            url = asset_directory + filename
            result = browser.evaluate('(async()=>{const response=await fetch(' + json.dumps(url) + ',{cache:"force-cache"});if(!response.ok)throw Error("Missing staged worker asset");const bytes=await response.arrayBuffer();const hash=await crypto.subtle.digest("SHA-256",bytes);return {length:bytes.byteLength,hash:Array.from(new Uint8Array(hash)).map(b=>b.toString(16).padStart(2,"0")).join("")}})()')
            relative = urlparse(url).path.lstrip('/')
            staged = QA_DIST / relative
            payload = staged.read_bytes()
            assert result == {'length': len(payload), 'hash': hashlib.sha256(payload).hexdigest()}, (filename, result)
        loader = browser.evaluate('(async()=>await(await fetch(' + json.dumps(asset_directory + NAMES[0]) + ',{cache:"force-cache"})).text())()')
        assert re.search(r'importScripts\(["\']\./analysis-worker\.js["\']\)', loader), loader
        assert re.search(r'wasm_bindgen\(["\']\./analysis-worker_bg\.wasm["\']\)', loader), loader
        # A missing/legacy meta URL must fail closed, never fall back to old URLs.
        before = browser.evaluate('window.__releaseWorkers.length')
        browser.evaluate('document.querySelector("meta[name=collagene-worker]").content="analysis-worker_loader.js"')
        mvp.click(browser, 'Run virtual digest')
        browser.wait('document.querySelector(".digest-status")?.textContent.includes("worker assets are missing")')
        assert browser.evaluate('window.__releaseWorkers.length') == before
        browser.evaluate('document.querySelector("meta[name=collagene-worker]").content=' + json.dumps(meta))
        inspect_catalog(browser)
        # Repeat with the SAME profile/cache. Startup errors remain accumulated.
        tape.events.clear()
        browser.call('Page.reload')
        browser.wait('document.body?.innerText.includes("Saved in this browser")')
        browser.wait('document.querySelector(".track-status")?.textContent.includes("Visible tracks")')
        inspect_catalog(browser)
        assert assert_worker_release(browser, tape, 'reload', isolated_tags=(19, 20)) == meta
        browser.screenshot('/tmp/collagene2-worker-release-desktop.png')
        print('PASS: existing HTTP cache/profile preserved; invalid old fixed URLs armed; full hashed loader/JS/WASM chain matches distribution; all constructors and actual child requests use one release directory; schema22 isolated autosave, retained tracks and catalog operations; missing meta fails closed; reload preserves working release; all startup errors checked')
    finally:
        browser.call('Page.removeScriptToEvaluateOnNewDocument', {'identifier': hook})
        browser.call('Fetch.disable')
        browser.call('Target.setAutoAttach', {'autoAttach': False, 'waitForDebuggerOnStart': False, 'flatten': True})
        browser.evaluate('if(window.__releaseOriginalWorker)window.Worker=window.__releaseOriginalWorker')


if __name__ == '__main__':
    main()
