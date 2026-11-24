#!/usr/bin/env python3
"""Real-source annotation lane/reflow acceptance on the existing QA browser.

Uses actual GenBank imports and worker reports. The 1 Mb / 10,000 annotation
fixture is bounded to a 120 s import budget; this is functional scale coverage,
not a machine-independent frame-rate or memory benchmark. Clears QA IndexedDB.
"""
import json
import time
from pathlib import Path

from browser_qa import Browser, QA_ORIGIN
from desktop_mvp_browser_qa import record, upload
from workspace_browser_qa import wait_file
from workspace_navigation_browser_qa import (
    ROWS, assert_unchanged, click_id, fill, persisted, rows, title,
)

DIRECTORY = Path('/tmp/collagene2-annotation-layout-qa')
CAP_FEATURE = 'select[aria-label="Annotation rows"]'
CAP_PRIMER = 'select[aria-label="Primer rows"]'
IMPORT_TIMEOUT = 120
SPARSE_BASES = 1_000_000
SPARSE_FEATURES = 10_000


def feature(label, location, kind='misc_feature'):
    return f'     {kind:<16}{location}\n                     /label="{label}"\n'


def button(b, label):
    b.pointer_click('Array.from(document.querySelectorAll("button")).find(e=>e.textContent.trim()===' + json.dumps(label) + ')')


def controls(b, opened):
    toggle = 'Array.from(document.querySelectorAll(".track-toolbar button")).find(e=>e.textContent.trim()==="Sequence tracks")'
    expanded = b.evaluate('(' + toggle + ').getAttribute("aria-expanded")')
    assert expanded in ('true', 'false'), ('Invalid expanded-state accessibility value', expanded)
    if expanded != str(opened).lower():
        b.pointer_click(toggle)
    b.wait(('!!' if opened else '!') + 'document.querySelector(' + json.dumps(CAP_FEATURE) + ')')
    assert b.evaluate('(' + toggle + ').getAttribute("aria-expanded")') == str(opened).lower()


def model(b):
    return b.evaluate('(()=>{const d=document.querySelector(".sequence-canvas").dataset;return {height:Number(d.rowHeight),features:Number(d.featureLanes),primers:Number(d.primerLanes),overflow:Number(d.layoutOverflow),source:d.layoutSource}})()')


def expect(b, height, features, primers, overflow, revision):
    expected = {'height': height, 'features': features, 'primers': primers, 'overflow': overflow, 'source': revision}
    b.wait('(()=>{const d=document.querySelector(".sequence-canvas")?.dataset;return d&&Number(d.rowHeight)===' + str(height) + '&&Number(d.featureLanes)===' + str(features) + '&&Number(d.primerLanes)===' + str(primers) + '&&Number(d.layoutOverflow)===' + str(overflow) + '&&d.layoutSource===' + json.dumps(revision) + '})()', timeout=60)
    assert model(b) == expected
    return expected


def ready(b, complete=True):
    text = "Visible tracks complete" if complete else "Visible tracks incomplete; inspect warnings"
    b.wait('document.querySelector(".track-status")?.textContent===' + json.dumps(text), timeout=60)


def scroll(b, top):
    b.evaluate('(()=>{const e=document.querySelector(".sequence-scroll");e.scrollTop=' + str(top) + ';e.dispatchEvent(new Event("scroll"));})()')
    b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')


def anchor(b):
    return b.evaluate('(()=>{const e=document.querySelector(".sequence-scroll"),c=document.querySelector(".sequence-canvas"),h=Number(c.dataset.rowHeight),cols=Math.max(10,Math.min(100,Math.floor((Math.max(200,e.clientWidth)-75)/8.5/10)*10));return {base:Math.floor(e.scrollTop/h)*cols,columns:cols,top:e.scrollTop}})()')


def viewport_report(b, complete=True):
    controls(b, True)
    disclosure = 'document.querySelector(".viewport-tracks")'
    if not b.evaluate('(' + disclosure + ').open'):
        b.pointer_click('(' + disclosure + ').querySelector("summary")')
    ready(b, complete=complete)
    path = DIRECTORY / 'viewport-tracks.json'
    path.unlink(missing_ok=True)
    button(b, 'Download viewport report')
    return json.loads(wait_file(path))


def pointer_canvas(b, x, y):
    # Canvas coordinates are CSS pixels. Verify a real hit before dispatching.
    b.evaluate('document.querySelector(".sequence-canvas").scrollIntoView({block:"nearest",inline:"nearest",behavior:"instant"})')
    point = b.evaluate('(()=>{const c=document.querySelector(".sequence-canvas"),r=c.getBoundingClientRect(),x=r.left+' + str(x) + ',y=r.top+' + str(y) + ';return {x,y,hit:document.elementFromPoint(x,y)===c}})()')
    assert point.pop('hit'), point
    b.call('Input.dispatchMouseEvent', {'type': 'mousePressed', 'button': 'left', 'clickCount': 1, **point})
    b.call('Input.dispatchMouseEvent', {'type': 'mouseReleased', 'button': 'left', 'clickCount': 1, **point})


def hold_messages(b):
    b.evaluate('(()=>{if(window.__layoutOriginalPost)throw Error("QA worker hold already installed");window.__layoutHeld=[];window.__layoutOriginalPost=Worker.prototype.postMessage;Worker.prototype.postMessage=function(...args){window.__layoutHeld.push([this,args]);};})()')


def release_messages(b):
    return b.evaluate('(()=>{const send=window.__layoutOriginalPost;if(!send)return 0;const pending=window.__layoutHeld.splice(0);Worker.prototype.postMessage=send;delete window.__layoutOriginalPost;for(const [worker,args]of pending)send.apply(worker,args);return pending.length})()')


def bottom_shrink(b, revision, document_id, rows_above_bottom):
    """Observe spacer shrink before any follow-up selection/scroll/retry action."""
    controls(b, True)
    fill(b, CAP_FEATURE, '8', 'change')
    expect(b, 192, 6, 1, 0, revision)
    b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')
    scroll(b, 1_000_000_000)
    old_max = b.evaluate('document.querySelector(".sequence-scroll").scrollTop')
    scroll(b, max(0, old_max - rows_above_bottom * 192))
    ready(b)
    position = anchor(b)
    # Choose an actual visible base line, avoiding all annotation hit boxes.
    row = int(position['top'] // 192)
    y = row * 192 - position['top'] + 30
    if y < 12:
        row += 1
        y += 192
    index = row * position['columns'] + 5
    assert 0 <= index < 4096 and y >= 12, (position, row, y, index)
    pointer_canvas(b, 64 + 5 * 8.5 + 4.25, y)
    literal = ('ACGT' * 1024)[index]
    selected_text = f'Selected bases {index+1} to {index+1}: {literal}'
    b.wait('document.querySelector(".sequence-scroll .sr-only")?.textContent===' + json.dumps(selected_text))
    before = anchor(b)
    before['height'] = model(b)['height']
    assert before['height'] == 192
    assert b.evaluate('document.querySelector("[role=tab][aria-selected=true]").dataset.documentId') == document_id

    fill(b, CAP_FEATURE, '1', 'change')
    expect(b, 102, 1, 1, 5, revision)
    # Read only from this point until geometry/selection assertions finish. In
    # particular, no ready/retry, scroll, click, or new selection can mask a
    # canceled anchor callback after the DOM clamps the old scroll position.
    b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(()=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))))')
    after = anchor(b)
    geometry = b.evaluate('(()=>{const e=document.querySelector(".sequence-scroll");return {maximum:e.scrollHeight-e.clientHeight,height:e.clientHeight}})()')
    fraction = (before['top'] % before['height']) / before['height']
    mapped = (before['base'] // after['columns'] + fraction) * 102
    expected = min(max(0, mapped), geometry['maximum'])
    assert abs(after['top'] - expected) <= 103, {'old': before, 'new': after, 'viewport': geometry, 'mapped': mapped, 'expected': expected, 'rows_above_bottom': rows_above_bottom}
    # The near-bottom case must exercise restoration below the new maximum;
    # otherwise a browser's automatic end clamp alone could satisfy the test.
    if rows_above_bottom:
        assert expected < geometry['maximum'] - 103, (before, after, geometry, expected)
    assert model(b) == {'height': 102, 'features': 1, 'primers': 1, 'overflow': 5, 'source': revision}
    assert b.evaluate('document.querySelector(".sequence-scroll .sr-only").textContent') == selected_text
    assert b.evaluate('document.querySelector("[role=tab][aria-selected=true]").dataset.documentId') == document_id
    print('PASS bottom shrink geometry: ' + json.dumps({'old_top': before['top'], 'new_top': after['top'], 'mapped': mapped, 'new_max': geometry['maximum'], 'expected': expected, 'rows_above_bottom': rows_above_bottom, 'selected_base': index+1}), flush=True)


def main():
    DIRECTORY.mkdir(exist_ok=True)
    small = 'ACGT' * 16
    fixtures = [
        ('layout_empty', small, False, ''),
        ('layout_source_only', small, False, feature('source_is_not_a_display_lane', '1..64', 'source')),
        ('layout_unresolved', small, False, feature('remote_location_retained', 'AB123456.1:5..9')),
        ('layout_one', small, False, feature('one', '2..12')),
        ('layout_distant', small, False, feature('left', '2..12') + feature('right', '40..52')),
        ('layout_join_gap', small, True, feature('wrapped', 'join(57..64,1..8)') + feature('in_gap', '20..30')),
        ('layout_overlap', 'ACGT' * 1024, False, ''.join(feature(f'overlap_{i}', '1..4096') for i in range(6)) + feature('primer_one', '2..8', 'primer_bind')),
        ('layout_primers', small, False, ''.join(feature(f'primer_{i}', '2..12', 'primer_bind') for i in range(3))),
    ]
    assert SPARSE_FEATURES * 100 == SPARSE_BASES
    sparse_annotations = ''.join(feature(f'sparse_{i:05d}', f'{i*100+1}..{i*100+5}') for i in range(SPARSE_FEATURES))
    fixtures.append(('layout_sparse_1Mb_10k', 'ACGT' * (SPARSE_BASES // 4), False, sparse_annotations))
    paths = {}
    for name, sequence, circular, annotations in fixtures:
        path = DIRECTORY / (name + '.gb')
        path.write_text(record(name, sequence, circular, annotations))
        assert path.stat().st_size < 8 * 1024 * 1024, 'QA fixture exceeded explicit file budget'
        paths[name] = path

    b = Browser()
    b.call('Browser.setDownloadBehavior', {'behavior': 'allow', 'downloadPath': str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride', {'width': 1600, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})
    b.call('Storage.clearDataForOrigin', {'origin': QA_ORIGIN, 'storageTypes': 'indexeddb'})
    b.exceptions.clear()
    b.console_errors.clear()  # Historical entries only, before the first fresh navigation.
    b.call('Page.navigate', {'url': QA_ORIGIN})
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    imports = {}
    for name, _, _, _ in fixtures:
        started = time.monotonic()
        upload(b, paths[name])
        b.wait('document.querySelector("h1")?.textContent===' + json.dumps(name), timeout=IMPORT_TIMEOUT)
        ready(b, complete=name != 'layout_unresolved')
        assert time.monotonic() - started < IMPORT_TIMEOUT, ('Import budget exceeded', name)
        imports[name] = next(row['id'] for row in rows(b) if row['name'] == name)
    baseline = persisted(b)
    revisions = {name: next(doc['head'] for doc in baseline[0]['documents'] if doc['id'] == document_id) for name, document_id in imports.items()}

    def activate(name):
        click_id(b, ROWS, imports[name])
        title(b, name)
        ready(b, complete=name != 'layout_unresolved')
        scroll(b, 0)

    for name, expected in [
        ('layout_empty', (66, 0, 0, 0)),
        ('layout_source_only', (66, 0, 0, 0)),
        ('layout_unresolved', (66, 0, 0, 0)),
        ('layout_one', (84, 1, 0, 0)),
        ('layout_distant', (84, 1, 0, 0)),
        ('layout_join_gap', (84, 1, 0, 0)),
    ]:
        activate(name)
        expect(b, *expected, revisions[name])

    activate('layout_unresolved')
    unresolved_report = viewport_report(b, complete=False)
    assert not unresolved_report['complete'] and not unresolved_report['items'], unresolved_report
    assert any("unresolved locations" in warning for warning in unresolved_report['warnings']), unresolved_report
    controls(b, False)
    b.click('Features')
    b.wait('document.querySelector(".inspector-content h2 .count")?.textContent==="1"')
    assert b.evaluate('document.querySelector(".feature-entry strong").textContent') == 'remote_location_retained'

    name = 'layout_overlap'
    activate(name)
    fixed = expect(b, 138, 3, 1, 3, revisions[name])
    assert b.evaluate('!!document.querySelector(".annotation-layout-note")')
    # Hidden lanes retain exact viewport evidence and native feature indices.
    report = viewport_report(b)
    annotations = [item for item in report['items'] if item['kind'] in ('Feature', 'Primer')]
    assert report['complete'] and len(annotations) == 7, report
    assert {item['source_feature'] for item in annotations} == set(range(7))
    assert {item['label'] for item in annotations} == {*(f'overlap_{i}' for i in range(6)), 'primer_one'}
    controls(b, False)
    button(b, 'Inspect annotations')
    b.wait('document.querySelector(".inspector-content h2")?.textContent.startsWith("Annotations")')
    assert b.evaluate('document.querySelectorAll(".feature-entry").length') == 7

    # Delayed real worker requests cannot determine the already-painted layout.
    # Queue and release unchanged outgoing messages; never fabricate a report.
    hold_messages(b)
    try:
        scroll(b, 6000)
        b.wait('window.__layoutHeld.length>0')
        assert model(b) == fixed
        b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')
        assert model(b) == fixed
    finally:
        released = release_messages(b)
    assert released > 0, 'No real worker request was delayed'
    ready(b)
    assert model(b) == fixed
    for top in (12000, 2000, 18000, 0):
        scroll(b, top)
        ready(b)
        assert model(b) == fixed

    # Reflow preserves the approximate first visible source base, independent of
    # the new pixel height. The source revision, IDs and native history stay fixed.
    scroll(b, 6000)
    before = anchor(b)
    controls(b, True)
    fill(b, CAP_FEATURE, '8', 'change')
    expect(b, 192, 6, 1, 0, revisions[name])
    b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')
    after = anchor(b)
    assert abs(after['base'] - before['base']) <= max(before['columns'], after['columns']), (before, after)
    assert not b.evaluate('!!document.querySelector(".annotation-layout-note")')
    fill(b, CAP_FEATURE, '1', 'change')
    expect(b, 102, 1, 1, 5, revisions[name])
    fill(b, CAP_FEATURE, '8', 'change')
    expect(b, 192, 6, 1, 0, revisions[name])
    for rows_above_bottom in (0, 8):
        bottom_shrink(b, revisions[name], imports[name], rows_above_bottom)
        assert_unchanged(baseline, persisted(b))
    fill(b, CAP_FEATURE, '8', 'change')
    expect(b, 192, 6, 1, 0, revisions[name])
    controls(b, False)
    scroll(b, 0)
    ready(b)
    # First actual annotation lane is at CSS y=66; select a painted bar.
    pointer_canvas(b, 64 + 4.25, 66)
    b.wait('document.querySelector(".track-selection strong")?.textContent==="overlap_0"')
    b.wait('document.querySelector(".annotation-details h3")?.textContent==="overlap_0"')
    assert b.evaluate('document.querySelector(".annotation-sequence").textContent') == 'ACGT' * 1024
    # Literal base6 is C. A real base hit after reflow must select exactly that
    # base, rather than an annotation lane or the old row's coordinate.
    pointer_canvas(b, 64 + 5 * 8.5 + 4.25, 30)
    b.wait('document.querySelector(".sequence-scroll .sr-only")?.textContent==="Selected bases 6 to 6: C"')
    controls(b, True)
    fill(b, CAP_FEATURE, '3', 'change')
    controls(b, False)

    activate('layout_primers')
    expect(b, 102, 0, 2, 1, revisions['layout_primers'])
    controls(b, True)
    fill(b, CAP_PRIMER, '4', 'change')
    expect(b, 120, 0, 3, 0, revisions['layout_primers'])
    fill(b, CAP_PRIMER, '1', 'change')
    expect(b, 84, 0, 1, 2, revisions['layout_primers'])
    fill(b, CAP_PRIMER, '2', 'change')
    controls(b, False)

    name = 'layout_sparse_1Mb_10k'
    activate(name)
    fixed = expect(b, 84, 1, 0, 0, revisions[name])
    b.click('Features')
    b.wait('document.querySelector(".inspector-content h2 .count")?.textContent==="10000"')
    # Verify actual million-base extent, not merely a small set of long labels.
    assert b.evaluate('document.querySelector(".document-meta").textContent.includes("1000000 bp")')
    for top in (50000, 250000, 1000000, 2000000, 0):
        scroll(b, top)
        ready(b)
        assert model(b) == fixed
    b.screenshot(DIRECTORY / 'sparse-1Mb-10k.png')
    assert_unchanged(baseline, persisted(b))
    # A fresh reload restores source evidence, not an injected display plan.
    b.call('Page.reload')
    title(b, name)
    ready(b)
    expect(b, 84, 1, 0, 0, revisions[name])
    assert_unchanged(baseline, persisted(b))
    assert not b.exceptions and not b.console_errors, (b.exceptions, b.console_errors)
    print('PASS: real GenBank whole-source lane plans; empty/source-only/distant/circular-gap packing; overflow retains all7 viewport annotations; feature/primer caps; delayed-worker/scroll stability; middle/bottom/near-bottom source-anchor and exact annotation/base selection after reflow; actual1Mb/10000 sparse annotations; immutable native source/history and reload; no fresh runtime errors', flush=True)


if __name__ == '__main__':
    main()
