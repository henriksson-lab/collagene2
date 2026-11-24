#!/usr/bin/env python3
"""Actual WASM run retention, exact source DAG and replay. Synthetic DNA only.
Run against the isolated Chrome :9224 / web :8080 used by browser_qa.py.
This clears this origin's IndexedDB; no application external tools are involved.
"""
import base64
import io
import json
import zipfile
from pathlib import Path

from browser_qa import QA_ORIGIN, Browser, mobile_enabled
from methods_browser_qa import record, fill, source, active, methods, compute, verify
from workspace_browser_qa import slots, wait_file, manifest

DIRECTORY = Path('/tmp/collagene2-run-history-qa')
ORIGIN = QA_ORIGIN


def close(b):
    b.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()')
    b.wait('!document.querySelector("[role=dialog]")')


def archive(b):
    p = DIRECTORY / 'collagene-workspace.cg2workspace'
    p.unlink(missing_ok=True)
    b.click('Workspace'); b.click('Save workspace file')
    data = wait_file(p); close(b)
    return data


def attachment_json(data, kind):
    """Read retained JSON independently of application decoding."""
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        result = []
        for entry in z.namelist():
            if not entry.startswith('blobs/'):
                continue
            try:
                value = json.loads(z.read(entry))
            except (ValueError, UnicodeDecodeError):
                continue
            if isinstance(value, dict) and 'run_id' in value and 'products' in value and 'plan' in value:
                result.append(value)
        return result


def record_dialog(b):
    b.click('Keep reconstructed products')
    b.wait('!!document.querySelector(".run-capture-panel")')


def details(b):
    b.click('History')
    b.wait('!!document.querySelector(".history-panel")')
    b.click('View cloning run')
    b.wait('!!document.querySelector(".run-details-panel")')


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for p in DIRECTORY.iterdir():
        if p.is_file(): p.unlink()
    b = Browser()
    b.call('Browser.setDownloadBehavior', {'behavior': 'allow', 'downloadPath': str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride', {'width': 1600, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})
    b.call('Storage.clearDataForOrigin', {'origin': ORIGIN, 'storageTypes': 'indexeddb'})
    b.call('Page.navigate', {'url': ORIGIN})
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    fixtures = [('ta_vector', 'ATGATATCGG', True), ('ta_vector', 'ATGATATCGG', True), ('pcr_template', 'ACGACCCCGGGGTAAC', False), ('ta_target', 'ATCGGATGATTACGACCCCGGGGTAACA', True)]
    b.click('Paste sequence'); fill(b, 'Sequence input', ''.join(record(*f) for f in fixtures)); b.click('Open sequence')
    b.wait('document.querySelector("h1")?.textContent==="ta_vector"')
    active(b, 'ta_target'); methods(b)
    source(b, 'Vector tailing source', 'ta_vector'); source(b, 'Insert tailing source', 'pcr_template')
    fill(b, 'Insert tailing forward primer', 'ACGA'); fill(b, 'Insert tailing reverse primer', 'GTTA')
    fill(b, 'Insert tailing terminal treatment', 'phosphorylate', 'change')
    compute(b); verify(b, 'ConsistentUnderAssumptions')
    before_data = archive(b); before = manifest(before_data)
    record_dialog(b)
    # Identical molecular revisions from different document/event identities
    # require an explicit choice instead of first-hash selection.
    selections = b.evaluate('Array.from(document.querySelectorAll(".run-source-binding select")).map(e=>({value:e.value,options:e.options.length}))')
    assert any(not e['value'] and e['options'] >= 3 for e in selections), selections
    before_slots = slots(b)
    b.click('Record cloning run')
    b.wait('document.querySelector(".run-capture-panel")?.textContent.includes("Choose the source event")')
    assert slots(b) == before_slots
    b.evaluate('Array.from(document.querySelectorAll(".run-source-binding select")).filter(e=>!e.value).forEach(e=>{e.selectedIndex=1;e.dispatchEvent(new Event("change",{bubbles:true}));})')
    assert all(b.evaluate('Array.from(document.querySelectorAll(".run-source-binding select")).map(e=>e.value)'))
    # Force unmount while its isolated worker is waiting for its request.
    # The cancellation interval must terminate it and unlock the unchanged UI.
    b.evaluate('window.qaPostMessage=Worker.prototype.postMessage;Worker.prototype.postMessage=function(){};true')
    b.click('Record cloning run')
    close(b)
    b.evaluate('Worker.prototype.postMessage=window.qaPostMessage;delete window.qaPostMessage;true')
    b.wait('!Array.from(document.querySelectorAll("[role=status]")).some(e=>e.textContent==="Updating workspace…")')
    assert slots(b) == before_slots, 'Unmounted recording changed storage'
    record_dialog(b)
    b.evaluate('Array.from(document.querySelectorAll(".run-source-binding select")).filter(e=>!e.value).forEach(e=>{e.selectedIndex=1;e.dispatchEvent(new Event("change",{bubbles:true}));})')
    # Synchronous second-put failure must abort the queued recovery-slot put.
    b.evaluate("""window.qaOriginalPut=IDBObjectStore.prototype.put;
      IDBObjectStore.prototype.put=function(value,key){if(key==='current')throw new DOMException('Synthetic current put failure','QuotaExceededError');return window.qaOriginalPut.call(this,value,key)};true""")
    b.click('Record cloning run')
    b.wait('document.querySelector(".run-capture-panel")?.textContent.includes("Run was not recorded")', timeout=60)
    b.evaluate('IDBObjectStore.prototype.put=window.qaOriginalPut;delete window.qaOriginalPut;true')
    assert slots(b) == before_slots, 'A failed transaction changed storage slots'
    assert b.evaluate('document.querySelector("h1").textContent') == 'ta_target'
    b.click('Record cloning run')
    b.wait('!document.querySelector(".run-capture-panel")', timeout=60)
    b.wait('document.querySelector("h1")?.textContent==="ta_target"')
    after_data = archive(b); after = manifest(after_data)
    runs = attachment_json(after_data, 'cloning-run-v1'); assert len(runs) == 1, runs
    run = runs[0]; products = run['products']; assert len(products) == len(run['plan']['steps']) and len(products) >= 5
    old_docs = {d['id']: d for d in before['documents']}
    new_docs = {d['id']: d for d in after['documents']}
    for doc_id, doc in old_docs.items(): assert new_docs[doc_id] == doc, (doc, new_docs[doc_id])
    assert after['active_document'] == before['active_document'] == run['target']['document_id']
    events = {e['id']: e for e in after['history']}
    inputs = {i['input_id']: i['event']['event_id'] for i in run['inputs']}
    for p, op in zip(products, run['plan']['steps']):
        expected = set()
        def sources(value):
            if isinstance(value, dict):
                if set(value) == {'Input'}: expected.add(inputs[value['Input']]); return
                if set(value) == {'Step'}: expected.add(products[value['Step']]['event_id']); return
                for item in value.values(): sources(item)
            elif isinstance(value, list):
                for item in value: sources(item)
        sources(op)
        assert set(events[p['event_id']]['parents']) == expected, (p, expected, events[p['event_id']])
        assert run['target']['event_id'] not in expected
    assert any(len(events[p['event_id']]['parents']) == 2 for p in products), 'Expected an actual two-source ligation join'
    recovery = slots(b)['previous']; assert recovery and manifest(base64.b64decode(recovery))['documents'] == before['documents']
    details(b)
    assert 'unverified on opening' in b.evaluate('document.querySelector(".run-details-panel").textContent')
    assert b.evaluate('document.querySelector(".run-replay-result")===null')
    b.click('Replay cloning run'); b.wait('!!document.querySelector(".run-replay-result")', timeout=60)
    assert 'products match recorded: true' in b.evaluate('document.querySelector(".run-replay-result").textContent')
    b.click('Download run replay JSON')
    replay = json.loads(wait_file(DIRECTORY / 'cloning-run-replay.json'))
    assert replay['products_match'] and replay['current_report']['status'] == 'ConsistentUnderAssumptions'
    close(b)
    b.call('Page.reload'); b.wait('document.body?.innerText.includes("Saved in this browser")')
    details(b); assert b.evaluate('document.querySelector(".run-replay-result")===null')
    close(b)
    # The portable archive restores historical records, never replay authority.
    b.click('Workspace')
    b.call('DOM.enable')
    node = b.call('DOM.getDocument')['root']['nodeId']
    file_node = b.call('DOM.querySelector', {'nodeId': node, 'selector': '[role=dialog] input[type=file][accept*=cg2workspace]'})['nodeId']
    b.call('DOM.setFileInputFiles', {'nodeId': file_node, 'files': [str(DIRECTORY / 'collagene-workspace.cg2workspace')]})
    b.wait('!document.querySelector("[role=dialog]")', timeout=60)
    details(b)
    assert b.evaluate('document.querySelector(".run-replay-result")===null')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride', {'width': 390, 'height': 844, 'deviceScaleFactor': 1, 'mobile': True})
        assert b.evaluate('document.documentElement.scrollWidth<=innerWidth'), 'Run details overflow mobile viewport'
    b.click('Replay cloning run'); b.wait('!!document.querySelector(".run-replay-result")', timeout=60)
    b.evaluate('document.querySelector(".modal").scrollTop=0')
    if mobile_enabled():
        b.screenshot('/tmp/collagene2-run-mobile.png')
    b.call('Emulation.setDeviceMetricsOverride', {'width': 1600, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})
    b.screenshot('/tmp/collagene2-run-desktop.png')
    assert not b.exceptions and not b.console_errors, (b.exceptions, b.console_errors)
    print('Cloning run actual product DAG, retained source heads, independent target, recovery, replay and reload/mobile checks passed'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__ == '__main__':
    main()
