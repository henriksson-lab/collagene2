#!/usr/bin/env python3
"""Real GenBank map labels, exact hit identity and source-preserving view modes.

Development-only acceptance on the already-running QA Chrome. Clears this origin's
IndexedDB. Actual 1 Mb / 10,000-feature input has explicit 120 s import and 60 s
rotation-stage budgets; this is not an FPS or memory-performance benchmark.
"""
import json
import math
import time
from pathlib import Path

from browser_qa import Browser, QA_ORIGIN
from desktop_mvp_browser_qa import record, upload
from workspace_navigation_browser_qa import ROWS, assert_unchanged, click_id, persisted, rows, title

DIRECTORY = Path('/tmp/collagene2-map-labels-qa')
EVIDENCE = '.map-stage canvas[data-label-total]'
MAX_HITS = 128  # Frozen renderer candidate/hit budget, not all source annotations.
IMPORT_TIMEOUT = 120
DENSE_COUNT = 10_000
DENSE_BASES = 1_000_000


def feature(label, location, kind='misc_feature'):
    return f'     {kind:<16}{location}\n                     /label="{label}"\n'


def settle(b):
    b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')


def visible(b, selector):
    return b.evaluate('(()=>{const e=document.querySelector(' + json.dumps(selector) + ');if(!e)return false;const r=e.getBoundingClientRect();return e.getClientRects().length>0&&r.width>1&&r.height>1&&getComputedStyle(e).visibility!=="hidden"})()')


def mode(b, label):
    selector = 'button[aria-label=' + json.dumps(label) + ']'
    b.pointer_click('document.querySelector(' + json.dumps(selector) + ')')
    b.wait('document.querySelector(' + json.dumps(selector) + ')?.getAttribute("aria-pressed")==="true"')
    settle(b)
    states = b.evaluate('Array.from(document.querySelectorAll("button[aria-label]" )).filter(e=>["Map view","Sequence view","Split view"].includes(e.getAttribute("aria-label"))).map(e=>({label:e.getAttribute("aria-label"),pressed:e.getAttribute("aria-pressed")}))')
    assert len(states) == 3 and sum(item['pressed'] == 'true' for item in states) == 1, states
    assert all(item['pressed'] in ('true', 'false') for item in states), states
    assert visible(b, '.map-panel') == (label != 'Sequence view'), (label, states)
    assert visible(b, '.sequence-panel') == (label != 'Map view'), (label, states)


def evidence(b, revision):
    b.wait('(()=>{const c=document.querySelector(' + json.dumps(EVIDENCE) + ');if(document.querySelectorAll(' + json.dumps(EVIDENCE) + ').length!==1||!c||c.dataset.labelSource!==' + json.dumps(revision) + ')return false;const r=c.getBoundingClientRect(),dpr=Math.min(devicePixelRatio,2);return r.width>1&&r.height>1&&Math.abs(c.width-r.width*dpr)<2&&Math.abs(c.height-r.height*dpr)<2&&!!c.dataset.labelHits})()', timeout=60)
    return b.evaluate('(()=>{const c=document.querySelector(' + json.dumps(EVIDENCE) + '),d=c.dataset,r=c.getBoundingClientRect();return {total:Number(d.labelTotal),displayed:Number(d.labelDisplayed),omitted:Number(d.labelOmitted),abbreviated:Number(d.labelAbbreviated),unplaced:Number(d.labelUnplaced),limited:d.labelLimited,source:d.labelSource,hits:JSON.parse(d.labelHits),width:r.width,height:r.height}})()')


def audit_rectangles(value, indices):
    for key in ('total', 'displayed', 'omitted', 'abbreviated', 'unplaced'):
        assert isinstance(value[key], int) and value[key] >= 0, (key, value)
    assert value['displayed'] + value['omitted'] == value['total'], value
    assert value['displayed'] == len(value['hits']) <= MAX_HITS, value
    assert value['abbreviated'] <= value['displayed'] and value['unplaced'] <= value['omitted'], value
    assert value['limited'] in ('true', 'false'), value
    if value['omitted'] or value['unplaced']:
        assert value['limited'] == 'true', value
    seen = set()
    for rect in value['hits']:
        assert isinstance(rect['index'], int) and rect['index'] in indices and rect['index'] not in seen, rect
        seen.add(rect['index'])
        assert all(isinstance(rect[key], (int, float)) and math.isfinite(rect[key]) for key in ('x', 'y', 'width', 'height')), rect
        assert rect['width'] > 0 and rect['height'] > 0, rect
        assert rect['x'] >= -.75 and rect['y'] >= -.75, (rect, value)
        assert rect['x'] + rect['width'] <= value['width'] + .75 and rect['y'] + rect['height'] <= value['height'] + .75, (rect, value)
    for i, first in enumerate(value['hits']):
        for second in value['hits'][i+1:]:
            dx = min(first['x']+first['width'], second['x']+second['width']) - max(first['x'], second['x'])
            dy = min(first['y']+first['height'], second['y']+second['height']) - max(first['y'], second['y'])
            assert dx <= .5 or dy <= .5, ('Overlapping visible label rectangles', first, second)


def click_label(b, index, revision):
    b.evaluate('document.querySelector(' + json.dumps(EVIDENCE) + ').scrollIntoView({block:"nearest",inline:"nearest",behavior:"instant"})')
    settle(b)
    value = evidence(b, revision)
    hit = next((item for item in value['hits'] if item['index'] == index), None)
    assert hit is not None, ('Required actual source feature has no label hit', index, value)
    point = b.evaluate('(()=>{const c=document.querySelector(' + json.dumps(EVIDENCE) + '),r=c.getBoundingClientRect(),x=r.left+' + str(hit['x']+hit['width']/2) + ',y=r.top+' + str(hit['y']+hit['height']/2) + ';return {x,y,unobscured:document.elementFromPoint(x,y)===c}})()')
    assert point.pop('unobscured'), point
    b.call('Input.dispatchMouseEvent', {'type': 'mousePressed', 'button': 'left', 'clickCount': 1, **point})
    b.call('Input.dispatchMouseEvent', {'type': 'mouseReleased', 'button': 'left', 'clickCount': 1, **point})


def rotate(b, angle):
    b.evaluate('document.querySelector(' + json.dumps(EVIDENCE) + ').scrollIntoView({block:"nearest",inline:"nearest",behavior:"instant"})')
    settle(b)
    circle = b.evaluate('(()=>{const c=document.querySelector(' + json.dumps(EVIDENCE) + '),r=c.getBoundingClientRect();return {x:r.left+r.width/2,y:r.top+r.height/2,r:Math.min(r.width,r.height)*.25}})()')
    start = {'x': circle['x'] + circle['r'], 'y': circle['y']}
    assert b.evaluate('document.elementFromPoint(' + str(start['x']) + ',' + str(start['y']) + ')===document.querySelector(' + json.dumps(EVIDENCE) + ')'), start
    b.call('Input.dispatchMouseEvent', {'type': 'mousePressed', 'button': 'left', 'buttons': 1, 'clickCount': 1, **start})
    point = start
    for step in range(1, 7):
        a = angle * step / 6
        point = {'x': circle['x'] + circle['r']*math.cos(a), 'y': circle['y'] + circle['r']*math.sin(a)}
        b.call('Input.dispatchMouseEvent', {'type': 'mouseMoved', 'button': 'left', 'buttons': 1, **point})
    b.call('Input.dispatchMouseEvent', {'type': 'mouseReleased', 'button': 'left', 'clickCount': 1, **point})
    settle(b)


def main():
    DIRECTORY.mkdir(exist_ok=True)
    sequence = 'ACGT' * 600
    spaced = ''.join(feature(f'f{i:02d}', f'{i*100+12}..{i*100+38}') for i in range(24))
    identity = feature('not_a_display_label', '1..2400', 'source') + feature('overlap_first', '100..140') + feature('overlap_second', '100..140') + feature('remote_unplaced', 'AB123456.1:5..9') + feature('separate', '1700..1740')
    long_name = 'Complete_annotation_name_' + 'long_label_' * 40
    oversized_name = 'Oversized_annotation_' + 'x' * 4100
    assert len(oversized_name.encode()) > 4096
    fixtures = [
        ('map_labels_circular', sequence, True, spaced),
        ('map_labels_linear', sequence, False, spaced),
        ('map_labels_identity', sequence, True, identity),
        ('map_labels_abbreviation', sequence, True, feature(long_name, '200..400')),
        ('map_labels_name_guard', sequence, True, feature('ordinary_retained', '100..140') + feature(oversized_name, '1000..1100')),
        ('map_labels_dense_10k', 'ACGT' * (DENSE_BASES//4), True, ''.join(feature(f'dense_{i:05d}', f'{i*100+1}..{i*100+12}') for i in range(DENSE_COUNT))),
    ]
    paths = {}
    for name, dna, circular, annotations in fixtures:
        path = DIRECTORY / (name + '.gb')
        path.write_text(record(name, dna, circular, annotations))
        assert path.stat().st_size < 8*1024*1024
        paths[name] = path

    b = Browser()
    b.call('Emulation.setDeviceMetricsOverride', {'width': 1600, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})
    b.call('Storage.clearDataForOrigin', {'origin': QA_ORIGIN, 'storageTypes': 'indexeddb'})
    b.exceptions.clear()
    b.console_errors.clear()  # Only historical entries, before fresh navigation.
    b.call('Page.navigate', {'url': QA_ORIGIN})
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    documents = {}
    for name, _, _, _ in fixtures:
        started = time.monotonic()
        upload(b, paths[name])
        b.wait('document.querySelector("h1")?.textContent===' + json.dumps(name), timeout=IMPORT_TIMEOUT)
        assert time.monotonic()-started < IMPORT_TIMEOUT, name
        documents[name] = next(row['id'] for row in rows(b) if row['name'] == name)
    baseline = persisted(b)
    revisions = {name: next(doc['head'] for doc in baseline[0]['documents'] if doc['id'] == document_id) for name, document_id in documents.items()}

    def activate(name):
        click_id(b, ROWS, documents[name])
        title(b, name)
        mode(b, 'Map view')
        return evidence(b, revisions[name])

    for name in ('map_labels_circular', 'map_labels_linear'):
        value = activate(name)
        audit_rectangles(value, set(range(24)))
        assert value['total'] == 24 and value['displayed'] > 12, value
        assert value['unplaced'] == 0
        b.screenshot(DIRECTORY / (name + '.png'))
        # Switching views changes presentation; the source/native record is fixed.
        for view in ('Sequence view', 'Split view', 'Map view'):
            mode(b, view)
            assert b.evaluate('document.querySelector("[role=tab][aria-selected=true]").dataset.documentId') == documents[name]
            title(b, name)
            if view != 'Sequence view':
                audit_rectangles(evidence(b, revisions[name]), set(range(24)))
        assert_unchanged(baseline, persisted(b))

    # Short f00 labels fit without variants. Install instrumentation only after
    # source installation, roomy sizing and a stable warm-up draw; it must not
    # turn expected initial measurements into a false rotation-cache failure.
    name = 'map_labels_circular'
    warm = activate(name)
    assert warm['displayed'] > 12 and warm['abbreviated'] == 0, warm
    rotate(b, .13)
    warm = evidence(b, revisions[name])
    assert warm['displayed'] > 12 and warm['abbreviated'] == 0, warm
    settle(b)
    b.evaluate('(()=>{if(window.__mapOriginalMeasureText)throw Error("QA measurement counter already installed");window.__mapMeasureCount=0;window.__mapOriginalMeasureText=CanvasRenderingContext2D.prototype.measureText;CanvasRenderingContext2D.prototype.measureText=function(...args){window.__mapMeasureCount++;return window.__mapOriginalMeasureText.apply(this,args);};})()')
    try:
        for angle in (.21, -.47, .83):
            rotate(b, angle)
            value = evidence(b, revisions[name])
            audit_rectangles(value, set(range(24)))
            assert (value['width'], value['height']) == (warm['width'], warm['height']), (warm, value)
            assert value['abbreviated'] == 0 and value['displayed'] > 12, value
            assert b.evaluate('window.__mapMeasureCount') == 0, 'Rotation remeasured unchanged short source labels'
    finally:
        b.evaluate('(()=>{const original=window.__mapOriginalMeasureText;if(original){CanvasRenderingContext2D.prototype.measureText=original;delete window.__mapOriginalMeasureText;}})()')
    assert_unchanged(baseline, persisted(b))

    # Source switching can happen while a pointer is held (for example from a
    # keyboard action). Do not send any incidental pointer event that would
    # cancel/release the old gesture before exercising the source-change guard.
    activate('map_labels_circular')
    b.click('Features')
    assert not b.evaluate('!!document.querySelector(".annotation-details,.feature-entry.active")')
    b.evaluate('document.querySelector(' + json.dumps(EVIDENCE) + ').scrollIntoView({block:"nearest",inline:"nearest",behavior:"instant"})')
    settle(b)
    pressed = b.evaluate('(()=>{const c=document.querySelector(' + json.dumps(EVIDENCE) + '),r=c.getBoundingClientRect(),x=r.left+r.width/2,y=r.top+r.height/2;return {x,y,hit:document.elementFromPoint(x,y)===c}})()')
    assert pressed.pop('hit'), pressed
    b.call('Input.dispatchMouseEvent', {'type': 'mousePressed', 'button': 'left', 'buttons': 1, 'clickCount': 1, **pressed})
    released = False
    try:
        target = 'map_labels_linear'
        # This explicit DOM click is essential to the regression: a normal row
        # pointer click would first release the gesture and hide the stale state.
        b.evaluate('Array.from(document.querySelectorAll(' + json.dumps(ROWS) + ')).find(e=>e.dataset.documentId===' + json.dumps(documents[target]) + ').click()')
        title(b, target)
        value = evidence(b, revisions[target])
        audit_rectangles(value, set(range(24)))
        assert value['hits'] and not b.evaluate('!!document.querySelector(".annotation-details,.feature-entry.active")')
        hit = value['hits'][0]
        release = b.evaluate('(()=>{const c=document.querySelector(' + json.dumps(EVIDENCE) + '),r=c.getBoundingClientRect(),x=r.left+' + str(hit['x']+hit['width']/2) + ',y=r.top+' + str(hit['y']+hit['height']/2) + ';return {x,y,hit:document.elementFromPoint(x,y)===c}})()')
        assert release.pop('hit'), release
        # No mouseMoved: stale drag=false would incorrectly select B's label.
        b.call('Input.dispatchMouseEvent', {'type': 'mouseReleased', 'button': 'left', 'clickCount': 1, **release})
        released = True
    finally:
        if not released:
            b.call('Input.dispatchMouseEvent', {'type': 'mouseReleased', 'button': 'left', 'clickCount': 1, **pressed})
    settle(b)
    assert not b.evaluate('!!document.querySelector(".annotation-details,.feature-entry.active")'), 'A source-A gesture selected a source-B annotation'
    assert b.evaluate('document.querySelector("[role=tab][aria-selected=true]").dataset.documentId') == documents[target]
    assert evidence(b, revisions[target])['source'] == revisions[target]
    assert_unchanged(baseline, persisted(b))

    name = 'map_labels_identity'
    value = activate(name)
    audit_rectangles(value, {1, 2, 4})
    assert value['total'] == 4 and value['unplaced'] == 1 and value['limited'] == 'true', value
    assert {1, 2} <= {hit['index'] for hit in value['hits']}, value
    # Identical coordinates have different source indices: selecting by span or
    # display ordinal would silently choose overlap_first instead of this label.
    click_label(b, 2, revisions[name])
    b.wait('document.querySelector(".annotation-details h3")?.textContent==="overlap_second"')
    assert b.evaluate('document.querySelector(".feature-entry.active strong")?.textContent') == 'overlap_second'
    assert b.evaluate('document.querySelector(".annotation-sequence").textContent') == sequence[99:140]
    audit_rectangles(evidence(b, revisions[name]), {1, 2, 4})

    name = 'map_labels_abbreviation'
    value = activate(name)
    audit_rectangles(value, {0})
    assert value['displayed'] == 1 and value['abbreviated'] == 1, value
    click_label(b, 0, revisions[name])
    b.wait('document.querySelector(".annotation-details h3")?.textContent===' + json.dumps(long_name))
    assert b.evaluate('document.querySelector(".annotation-sequence").textContent') == sequence[199:400]

    name = 'map_labels_name_guard'
    value = activate(name)
    audit_rectangles(value, {0, 1})
    assert value['total'] == 2 and value['displayed'] == 0 and value['omitted'] == 2, value
    assert value['limited'] == 'true' and value['abbreviated'] == 0 and value['unplaced'] == 0, value
    assert not value['hits'], 'Oversized source name must withhold the whole map label cache, not publish a prefix'
    b.click('Features')
    b.wait('document.querySelector(".inspector-content h2 .count")?.textContent==="2"')
    actual_names = b.evaluate('Array.from(document.querySelectorAll(".feature-entry strong")).map(e=>e.textContent)')
    assert actual_names == ['ordinary_retained', oversized_name], 'Full source names must survive map cache limits'
    assert b.evaluate('(()=>{const e=document.querySelector(".inspector");return e.scrollWidth<=e.clientWidth+1})()'), 'Long imported names must wrap inside the inspector'
    assert_unchanged(baseline, persisted(b))

    name = 'map_labels_dense_10k'
    value = activate(name)
    assert b.evaluate('document.querySelector(".document-meta").textContent.includes("1000000 bp")')
    b.click('Features')
    b.wait('document.querySelector(".inspector-content h2 .count")?.textContent==="10000"')
    audit_rectangles(value, set(range(DENSE_COUNT)))
    assert value['total'] == DENSE_COUNT and value['omitted'] > 0 and value['limited'] == 'true', value
    assert value['unplaced'] == 0 and value['displayed'] > 0
    initial_hits = value['hits']
    started = time.monotonic()
    for angle in (.21, -.47, .83, .37):
        rotate(b, angle)
        value = evidence(b, revisions[name])
        audit_rectangles(value, set(range(DENSE_COUNT)))
        assert value['total'] == DENSE_COUNT and value['limited'] == 'true' and value['omitted'] > 0, value
        assert time.monotonic()-started < 60, 'Bounded dense rotation-stage budget exceeded'
    assert value['hits'] != initial_hits, 'Rotation did not update label hit geometry'
    b.screenshot(DIRECTORY / 'dense-10k-rotated.png')
    assert_unchanged(baseline, persisted(b))
    b.call('Page.reload')
    title(b, name)
    mode(b, 'Map view')
    value = evidence(b, revisions[name])
    audit_rectangles(value, set(range(DENSE_COUNT)))
    assert value['total'] == DENSE_COUNT and value['omitted'] > 0
    assert_unchanged(baseline, persisted(b))
    assert not b.exceptions and not b.console_errors, (b.exceptions, b.console_errors)
    print('PASS: actual24-feature circular/linear maps exceed12 labels when roomy; bounded unclipped nonoverlapping text hits; source-switch held-gesture cancellation; same-span exact original-feature selection; full-name inspection after abbreviation and whole-cache oversized-name guard; zero rotation text remeasurement; actual1Mb/10000-feature honest omissions and bounded rotation; three view modes and immutable native source/history/reload; no fresh runtime errors', flush=True)


if __name__ == '__main__':
    main()
