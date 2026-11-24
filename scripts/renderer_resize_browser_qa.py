#!/usr/bin/env python3
"""Desktop canvas viewport regression against an already-running CDP browser.

Only fixture imports and QA-local CSS change the page. In particular, the
scrollbar test sends no synthetic resize/scroll/input event to repair the canvas.
"""
import json
import math
import time
from pathlib import Path

from browser_qa import Browser, QA_ORIGIN
from desktop_mvp_browser_qa import document, record, upload
from workspace_navigation_browser_qa import assert_unchanged, saved_snapshot


DIRECTORY = Path('/tmp/collagene2-renderer-resize-qa')
SHORT = 'GATTACACCGTTAAGCTGCATACGATTCGAGCTTAGGCTAACGTACCATGCTAGCTTACGGATC'
LONG = ('ACGTTGCAAGCTTAGGATCCATGCTACGATTCGGAATTCACCGTACGATTAGCGTACCATGCA' * 40)
METRICS = '''(()=>{
 const e=document.querySelector('.sequence-scroll'),c=e.querySelector('canvas');
 const r=e.getBoundingClientRect(),cr=c.getBoundingClientRect();
 return {borderWidth:r.width,borderHeight:r.height,clientWidth:e.clientWidth,
 clientHeight:e.clientHeight,offsetWidth:e.offsetWidth,
 usableWidth:r.width-(e.offsetWidth-e.clientWidth),usableHeight:r.height-(e.offsetHeight-e.clientHeight),scrollWidth:e.scrollWidth,
 scrollHeight:e.scrollHeight,scrollTop:e.scrollTop,canvasWidth:cr.width,
 canvasHeight:cr.height,backingWidth:c.width,backingHeight:c.height,
 dpr:Math.min(devicePixelRatio,2),documentWidth:document.documentElement.scrollWidth,
 viewport:innerWidth,active:document.querySelector('[role=tab][aria-selected=true]')?.dataset.documentId};
})()'''


def frames(b, count=5):
    b.evaluate('new Promise(resolve=>{let left=' + str(count) + ';function step(){if(--left)requestAnimationFrame(step);else resolve();}requestAnimationFrame(step);})')


def correct(m):
    return (abs(m['canvasWidth'] - m['usableWidth']) <= .1
            and abs(m['canvasHeight'] - m['usableHeight']) <= .1
            and m['backingWidth'] == math.ceil(m['usableWidth'] * m['dpr'])
            and m['backingHeight'] == math.ceil(m['usableHeight'] * m['dpr'])
            and m['scrollWidth'] <= m['clientWidth'] + 1
            and m['documentWidth'] <= m['viewport'] + 1)


def settled(b, label, observations):
    # Reading metrics / advancing RAF never dispatches an application event.
    frames(b)
    deadline = time.monotonic() + 8
    while True:
        m = b.evaluate(METRICS)
        if correct(m):
            observations.append({'case': label, **m})
            return m
        if time.monotonic() >= deadline:
            observations.append({'case': label, **m})
            raise AssertionError((label, m))
        frames(b, 2)
        time.sleep(.05)


def viewport(b, width):
    b.call('Emulation.setDeviceMetricsOverride', {
        'width': width, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})


def selected_literal(b, sequence, document_id):
    m = b.evaluate(METRICS)
    assert m['active'] == document_id and m['scrollTop'] == 0, m
    columns = min(100, max(10, int((m['clientWidth'] - 75) // 8.5) // 10 * 10))
    index = min(columns, len(sequence)) - 1
    point = b.evaluate('(()=>{const c=document.querySelector(".sequence-canvas"),r=c.getBoundingClientRect();return {x:r.x+68+' + str(index) + '*8.5,y:r.y+25};})()')
    assert b.evaluate('document.elementFromPoint(' + str(point['x']) + ',' + str(point['y']) + ')===document.querySelector(".sequence-canvas")'), point
    for kind in ('mousePressed', 'mouseReleased'):
        b.call('Input.dispatchMouseEvent', {'type': kind, 'button': 'left', 'clickCount': 1, **point})
    expected = f'Selected bases {index + 1} to {index + 1}: {sequence[index]}'
    b.wait('Array.from(document.querySelectorAll(".sequence-scroll .sr-only")).some(e=>e.textContent===' + json.dumps(expected) + ')')
    assert b.evaluate(METRICS)['active'] == document_id


def main():
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, sequence in [('RESIZE_SHORT', SHORT), ('RESIZE_LONG', LONG)]:
        path = DIRECTORY / (name + '.gb')
        path.write_text(record(name, sequence, False))
        paths[name] = path
    b = Browser()
    observations = []
    try:
        viewport(b, 1600)
        b.call('Storage.clearDataForOrigin', {'origin': QA_ORIGIN, 'storageTypes': 'indexeddb'})
        # Clear historical Runtime.enable replay before navigation, never after startup.
        b.exceptions.clear()
        b.console_errors.clear()
        b.call('Page.navigate', {'url': QA_ORIGIN})
        b.wait('document.body?.innerText.includes("Saved in this browser")')
        for name in ('RESIZE_SHORT', 'RESIZE_LONG'):
            upload(b, paths[name])
            b.wait('document.querySelector("h1")?.textContent===' + json.dumps(name))
            b.wait('document.querySelector(".track-status")?.textContent==="Visible tracks complete"')
        before = saved_snapshot(b)
        for path in paths.values():
            assert any(path.read_bytes() in blob for blob in before[1].values()), ('Original fixture evidence missing', path.name)
        source_id = b.evaluate(METRICS)['active']
        # Finish source-bound worker rendering before the no-event probe. Otherwise
        # a late scientific response could redraw and conceal the stale-size bug.
        time.sleep(.4)
        frames(b)
        assert b.evaluate(METRICS)['scrollHeight'] > b.evaluate(METRICS)['clientHeight']
        b.evaluate('document.querySelector(".sequence-scroll").style.overflowY="hidden"')
        # Establish the hidden-scrollbar baseline using a real viewport resize.
        viewport(b, 1580)
        settled(b, 'hidden baseline resized', observations)
        viewport(b, 1600)
        baseline = settled(b, 'hidden baseline', observations)
        for mode in ('scroll', 'hidden'):
            b.evaluate('document.querySelector(".sequence-scroll").style.overflowY=' + json.dumps(mode))
            frames(b)
            current = b.evaluate(METRICS)
            assert abs(current['borderWidth'] - baseline['borderWidth']) < .1, (baseline, current)
            assert abs(current['borderHeight'] - baseline['borderHeight']) < .1, (baseline, current)
            assert current['active'] == source_id and current['scrollTop'] == 0, current
            if mode == 'scroll' and current['clientWidth'] == baseline['clientWidth']:
                print('SKIP scrollbar-only width delta: platform uses overlay scrollbars; all other resize checks remain required', flush=True)
            settled(b, 'no-event overflow-y ' + mode, observations)
        b.evaluate('document.querySelector(".sequence-scroll").style.removeProperty("overflow-y")')
        settled(b, 'restored natural scrollbar', observations)
        for width in (1366, 1200, 1600):
            viewport(b, width)
            settled(b, 'viewport ' + str(width), observations)
            for state in ('collapsed', 'expanded'):
                b.pointer_click('document.querySelector(".workspace-rail button[aria-label=Documents]")')
                settled(b, 'rail ' + state + ' at ' + str(width), observations)
            selected_literal(b, LONG, source_id)
        # Natural document-length changes, with original identities retained.
        for name, sequence in [('RESIZE_SHORT', SHORT), ('RESIZE_LONG', LONG)]:
            document(b, name)
            b.wait('document.querySelector(".track-status")?.textContent==="Visible tracks complete"')
            m = settled(b, 'natural switch ' + name, observations)
            if name == 'RESIZE_SHORT':
                assert m['offsetWidth'] == m['clientWidth'] and abs(m['borderHeight']-m['usableHeight']) < .1, ('Short record retained unnecessary scrollbars',m)
                assert m['scrollHeight'] <= m['clientHeight'] + 1, m
                # Grow only a fraction of a CSS pixel: rounded client sizes stay
                # unchanged, yet the actual canvas box must follow the observer.
                b.evaluate('document.querySelector(".sequence-scroll").style.width=' + json.dumps(str(m['clientWidth'] - .25) + 'px'))
                baseline=settled(b,'fractional viewport baseline',observations)
                b.evaluate('document.querySelector(".sequence-scroll").style.width=' + json.dumps(str(baseline['borderWidth'] + .125) + 'px'))
                frames(b)
                probe=b.evaluate(METRICS)
                assert probe['clientWidth']==baseline['clientWidth'] and probe['clientHeight']==baseline['clientHeight'], (baseline,probe)
                assert probe['borderWidth']>baseline['borderWidth'], (baseline,probe)
                settled(b,'fractional-only viewport growth',observations)
                b.evaluate('document.querySelector(".sequence-scroll").style.removeProperty("width")')
                m=settled(b,'fractional viewport restored',observations)
            selected_literal(b, sequence, m['active'])
        assert_unchanged(before, saved_snapshot(b))
        frames(b)
        assert not b.exceptions, b.exceptions
        assert not b.console_errors, b.console_errors
        b.screenshot(DIRECTORY / 'renderer-resize-desktop.png')
        print('PASS: no-event scrollbar-only canvas resizing; desktop rail/window widths; literal last-visible-column base identity; short/long imports and unchanged source/history; no startup/runtime errors', flush=True)
    finally:
        (DIRECTORY / 'metrics.json').write_text(json.dumps(observations, indent=2))
        b.socket.close()


if __name__ == '__main__':
    main()
