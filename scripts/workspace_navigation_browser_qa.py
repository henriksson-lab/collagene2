#!/usr/bin/env python3
"""Desktop document navigation acceptance; requires the existing isolated CDP browser.

Clears only QA_ORIGIN's IndexedDB. Uses actual imports, pointer/keyboard events and
saved native workspace evidence; no application state is installed by the harness.
"""
import base64
import io
import json
import time
import zipfile
from pathlib import Path

from browser_qa import QA_ORIGIN, Browser
from desktop_mvp_browser_qa import record, upload
from workspace_browser_qa import slots

DIRECTORY = Path('/tmp/collagene2-workspace-navigation-qa')
FIXTURES = [
    ('Zulu_plasmid', 'ACGT' * 31, True),
    ('alpha_vector', 'GATTACA' * 9, False),
    ('Middle_sequence_with_a_deliberately_long_document_identifier', 'AACCGGTT' * 23, True),
]


def fill(b, selector, value, event='input'):
    b.pointer_click('document.querySelector(' + json.dumps(selector) + ')')
    b.evaluate('(()=>{const e=document.querySelector(' + json.dumps(selector) + ');e.value=' + json.dumps(value) + ';e.dispatchEvent(new Event(' + json.dumps(event) + ',{bubbles:true}));})()')


def press(b, key, code=None):
    b.call('Input.dispatchKeyEvent', {'type': 'keyDown', 'key': key, 'code': code or key})
    b.call('Input.dispatchKeyEvent', {'type': 'keyUp', 'key': key, 'code': code or key})


def current(b):
    return b.evaluate('document.querySelector("h1")?.textContent')


def title(b, name):
    b.wait('document.querySelector("h1")?.textContent===' + json.dumps(name))


def saved_snapshot(b):
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    encoded = slots(b)['current']
    assert encoded, 'No persisted native workspace'
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(encoded))) as archive:
        manifest = json.loads(archive.read('manifest.json'))
        blobs = {name: archive.read(name) for name in archive.namelist() if name.startswith('blobs/')}
    # Active-document changes are expected; revisions/history and original source
    # payloads must remain identical throughout navigation and reload.
    return manifest, blobs


def assert_unchanged(before, after):
    old, old_blobs = before
    new, new_blobs = after
    for key in ('documents', 'history', 'revisions', 'blobs', 'candidates', 'attachments'):
        assert old[key] == new[key], ('Navigation changed native evidence', key)
    assert old_blobs == new_blobs, 'Navigation changed source/revision blobs'


def viewport(b, width, height=1000):
    b.call('Emulation.setDeviceMetricsOverride', {'width': width, 'height': height, 'deviceScaleFactor': 1, 'mobile': False})
    b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')
    assert b.evaluate('document.documentElement.scrollWidth<=innerWidth+1 && document.body.scrollWidth<=innerWidth+1'), b.evaluate('({width:innerWidth,html:document.documentElement.scrollWidth,body:document.body.scrollWidth})')


EXPLORER = 'aside[aria-label="Document explorer"]'
TABS = '[role="tablist"][aria-label="Open sequence documents"]'
SEARCH = 'input[aria-label="Filter documents"]'
SORT = 'select[aria-label="Sort documents"]'
ROWS = EXPLORER + ' .document-entry[data-document-id]'
TAB_BUTTONS = TABS + ' [role="tab"][data-document-id]'


def rows(b):
    return b.evaluate('Array.from(document.querySelectorAll(' + json.dumps(ROWS) + ')).map(e=>({id:e.dataset.documentId,name:e.querySelector("strong").textContent,metadata:e.querySelector("small").textContent}))')


def selected_id(b):
    return b.evaluate('document.querySelector(' + json.dumps(TAB_BUTTONS + '[aria-selected="true"]') + ')?.dataset.documentId')


def selected(b, document_id, name, focused=False):
    title(b, name)
    b.wait('document.querySelector(' + json.dumps(TAB_BUTTONS + '[aria-selected="true"]') + ')?.dataset.documentId===' + json.dumps(document_id))
    assert b.evaluate('document.querySelectorAll(' + json.dumps(TAB_BUTTONS + '[aria-selected="true"]') + ').length===1')
    assert b.evaluate('document.querySelectorAll(' + json.dumps(TAB_BUTTONS + '[tabindex="0"]') + ').length===1')
    assert b.evaluate('document.querySelector(' + json.dumps(ROWS + '[aria-current="page"]') + ')?.dataset.documentId===' + json.dumps(document_id))
    if focused:
        b.wait('document.activeElement?.dataset.documentId===' + json.dumps(document_id) + ' && document.activeElement?.getAttribute("role")==="tab"')


def click_id(b, selector, document_id):
    b.pointer_click('Array.from(document.querySelectorAll(' + json.dumps(selector) + ')).find(e=>e.dataset.documentId===' + json.dumps(document_id) + ')')


def persisted(b):
    # Wait for the native active ID too: an old "Saved" label is insufficient
    # while the next navigation transaction is queued.
    deadline = time.monotonic() + 30
    while time.monotonic() < deadline:
        snapshot = saved_snapshot(b)
        if snapshot[0]['active_document'] == selected_id(b):
            return snapshot
        time.sleep(.1)
    raise AssertionError('Selected document did not persist atomically')



def pagination(b, before):
    count_before = len(before[0]['documents'])
    path = DIRECTORY / 'paged-documents.gb'
    path.write_text(''.join(record(f'page_{i:03}', 'ACGTACGT', False) for i in range(101)))
    upload(b, path)
    title(b, 'page_000')
    total = count_before + 101
    b.wait('document.querySelectorAll(' + json.dumps(TAB_BUTTONS) + ').length===' + str(total))
    b.wait('document.querySelectorAll(' + json.dumps(ROWS) + ').length===100')
    first_id = selected_id(b)
    click_id(b, TAB_BUTTONS, first_id)
    press(b, 'End')
    title(b, 'page_100')
    last_id = selected_id(b)
    selected(b, last_id, 'page_100', focused=True)
    assert len(rows(b)) == total - 100
    assert b.evaluate('document.querySelector(' + json.dumps(TABS) + ').scrollLeft>0'), 'Active end tab was not revealed'
    b.pointer_click('document.querySelector(' + json.dumps('button[aria-label="Previous documents"]') + ')')
    b.wait('document.querySelectorAll(' + json.dumps(ROWS) + ').length===100')
    assert selected_id(b) == last_id and current(b) == 'page_100', 'Paging changed active document'
    middle = next(row for row in rows(b) if row['name'] == 'page_010')
    click_id(b, ROWS, middle['id'])
    selected(b, middle['id'], middle['name'])
    b.pointer_click('document.querySelector(' + json.dumps('button[aria-label="Next documents"]') + ')')
    b.wait('document.querySelectorAll(' + json.dumps(ROWS) + ').length===' + str(total - 100))
    assert selected_id(b) == middle['id'] and current(b) == middle['name']
    click_id(b, ROWS, last_id)
    selected(b, last_id, 'page_100')
    after = persisted(b)
    for doc in before[0]['documents']:
        assert doc in after[0]['documents'], 'Paging/import changed existing document heads/history'
    for event in before[0]['history']:
        assert event in after[0]['history'], 'Paging/import discarded existing history'
    assert all(after[1].get(name) == blob for name, blob in before[1].items()), 'Paging/import discarded historical source evidence'
    b.call('Page.reload')
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    selected(b, last_id, 'page_100')
    assert len(rows(b)) == total - 100
    assert_unchanged(after, persisted(b))


def main():
    DIRECTORY.mkdir(exist_ok=True)
    paths = []
    for name, sequence, circular in FIXTURES:
        path = DIRECTORY / (name + '.gb')
        path.write_text(record(name, sequence, circular, '     misc_feature    2..8\n                     /label="Navigation source marker"\n'))
        paths.append(path)
    b = Browser()
    b.call('Emulation.setDeviceMetricsOverride', {'width': 1600, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})
    b.call('Storage.clearDataForOrigin', {'origin': QA_ORIGIN, 'storageTypes': 'indexeddb'})
    # Discard only historical errors delivered by Runtime.enable, before fresh navigation.
    b.exceptions.clear()
    b.console_errors.clear()
    b.call('Page.navigate', {'url': QA_ORIGIN})
    b.wait('!!document.querySelector(' + json.dumps(EXPLORER) + ')')
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    initial_count = len(rows(b))
    for (name, _, _), path in zip(FIXTURES, paths):
        b.pointer_click('document.querySelector(' + json.dumps(EXPLORER + ' button[aria-label="Open sequence file"]') + ')')
        upload(b, path)
        title(b, name)
    b.wait('document.querySelectorAll(' + json.dumps(ROWS) + ').length===' + str(initial_count + len(FIXTURES)))
    original = rows(b)
    by_name = {row['name']: row['id'] for row in original}
    assert len(by_name) == len(original), original
    for name, sequence, circular in FIXTURES:
        row = next(row for row in original if row['name'] == name)
        assert row['metadata'].startswith(f'{len(sequence)} bp · {"Circular" if circular else "Linear"} DNA'), row
    baseline = persisted(b)
    for path in paths:
        assert any(path.read_bytes() in blob for blob in baseline[1].values()), ('Original GenBank missing from native source evidence', path.name)

    # Search both by case-insensitive name and stable document identity. An empty
    # result must leave the independently selected editor/tab unchanged.
    fill(b, SEARCH, 'ALPHA_VECTOR')
    b.wait('document.querySelectorAll(' + json.dumps(ROWS) + ').length===1')
    assert rows(b)[0]['id'] == by_name['alpha_vector']
    click_id(b, ROWS, by_name['alpha_vector'])
    selected(b, by_name['alpha_vector'], 'alpha_vector')
    fill(b, SEARCH, 'no-such-navigation-document')
    b.wait('document.querySelector(".explorer-empty")?.textContent==="No matching documents."')
    assert current(b) == 'alpha_vector' and selected_id(b) == by_name['alpha_vector']
    fill(b, SEARCH, by_name['Zulu_plasmid'])
    b.wait('document.querySelectorAll(' + json.dumps(ROWS) + ').length===1')
    assert rows(b)[0]['name'] == 'Zulu_plasmid'
    fill(b, SEARCH, '')
    b.wait('document.querySelectorAll(' + json.dumps(ROWS) + ').length===' + str(len(original)))
    fill(b, SORT, 'name', 'change')
    b.wait('Array.from(document.querySelectorAll(' + json.dumps(ROWS + ' strong') + ')).map(e=>e.textContent).join("|")===' + json.dumps('|'.join(sorted(by_name, key=str.lower))))
    name_order = [row['id'] for row in rows(b)]
    fill(b, SORT, 'length', 'change')
    b.wait('Array.from(document.querySelectorAll(' + json.dumps(ROWS + ' small') + ')).map(e=>Number(e.textContent.split(" ")[0])).every((n,i,a)=>!i||a[i-1]>=n)')
    fill(b, SORT, 'source', 'change')
    # These local fixtures and the demo have no provider declaration: Source
    # sorting must use the documented name fallback rather than fabricate one.
    b.wait('Array.from(document.querySelectorAll(' + json.dumps(ROWS) + ')).map(e=>e.dataset.documentId).join("|")===' + json.dumps('|'.join(name_order)))
    fill(b, SORT, 'workspace', 'change')
    b.wait('Array.from(document.querySelectorAll(' + json.dumps(ROWS) + ')).map(e=>e.dataset.documentId).join("|")===' + json.dumps('|'.join(row['id'] for row in original)))

    for row in original:
        click_id(b, ROWS, row['id'])
        selected(b, row['id'], row['name'])
    # Tab order stays workspace order regardless of explorer sorting/filtering.
    tabs = b.evaluate('Array.from(document.querySelectorAll(' + json.dumps(TAB_BUTTONS) + ')).map(e=>({id:e.dataset.documentId,name:e.querySelector("span").textContent}))')
    assert [tab['id'] for tab in tabs] == [row['id'] for row in original]
    click_id(b, TAB_BUTTONS, tabs[1]['id'])
    selected(b, tabs[1]['id'], tabs[1]['name'], focused=True)
    for key, index in [('ArrowRight', 2), ('Home', 0), ('ArrowLeft', len(tabs)-1), ('ArrowRight', 0), ('End', len(tabs)-1), ('ArrowLeft', len(tabs)-2)]:
        press(b, key)
        selected(b, tabs[index]['id'], tabs[index]['name'], focused=True)
    chosen = tabs[-2]
    assert_unchanged(baseline, persisted(b))

    toggle = 'document.querySelector(' + json.dumps('.workspace-rail button[aria-label="Documents"]') + ')'
    for width in (1600, 1280, 1200):
        viewport(b, width)
        assert b.evaluate('document.querySelector(' + json.dumps(EXPLORER) + ').getBoundingClientRect().width<=280'), 'Explorer is not compact'
        b.pointer_click(toggle)
        b.wait(toggle + '.getAttribute("aria-expanded")==="false"')
        assert b.evaluate('getComputedStyle(document.querySelector(' + json.dumps(EXPLORER) + ')).display==="none"')
        assert selected_id(b) == chosen['id'] and current(b) == chosen['name']
        viewport(b, width)
        b.pointer_click(toggle)
        b.wait(toggle + '.getAttribute("aria-expanded")==="true"')
        selected(b, chosen['id'], chosen['name'])
        viewport(b, width)
        b.screenshot(DIRECTORY / f'navigation-{width}.png')
    # A laptop's vertical space matters as much as width. Keep the document
    # status and working map visible without scrolling the whole application.
    for width, height in ((1366, 768), (1280, 800)):
        viewport(b, width, height)
        selected(b, chosen['id'], chosen['name'])
        geometry = b.evaluate('(()=>{const footer=document.querySelector(".statusbar").getBoundingClientRect(),map=document.querySelector(".map-stage").getBoundingClientRect(),list=document.querySelector(".explorer-document-list").getBoundingClientRect();return {footerTop:footer.top,footerBottom:footer.bottom,mapWidth:map.width,mapHeight:map.height,listBottom:list.bottom,height:innerHeight};})()')
        assert 0 <= geometry['footerTop'] < geometry['footerBottom'] <= height + 1, geometry
        assert geometry['mapWidth'] > 0 and geometry['mapHeight'] >= 200, geometry
        assert geometry['listBottom'] <= height + 1, geometry
        b.screenshot(DIRECTORY / f'navigation-{width}x{height}.png')
    before_reload = persisted(b)
    assert_unchanged(baseline, before_reload)
    b.call('Page.reload')
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    selected(b, chosen['id'], chosen['name'])
    assert_unchanged(baseline, persisted(b))
    assert selected_id(b) == before_reload[0]['active_document']
    viewport(b, 1200)

    # Intentional edit follows the immutable navigation phase. Stable keyed rows
    # and tabs must update their labels/metadata without replacing document IDs.
    before_edit = persisted(b)
    previous_doc = next(doc for doc in before_edit[0]['documents'] if doc['id'] == chosen['id'])
    row_expression = 'Array.from(document.querySelectorAll(' + json.dumps(ROWS) + ')).find(e=>e.dataset.documentId===' + json.dumps(chosen['id']) + ')'
    tab_expression = 'Array.from(document.querySelectorAll(' + json.dumps(TAB_BUTTONS) + ')).find(e=>e.dataset.documentId===' + json.dumps(chosen['id']) + ')'
    old_metadata = b.evaluate('(' + row_expression + ').querySelector("small").textContent')
    old_length = b.evaluate('(' + row_expression + ').dataset.documentLength')
    renamed = chosen['name'] + '_renamed'
    # Use the file menu's visible disclosure, then real pointer targets.
    more = 'document.querySelector(' + json.dumps('summary[aria-label="More file actions"]') + ')'
    if not b.evaluate('(' + more + ').parentElement.open'):
        b.pointer_click(more)
    b.pointer_click('Array.from(document.querySelectorAll(".mvp-file-menu button")).find(e=>e.textContent==="Document settings")')
    fill(b, 'input[aria-label="Document name"]', renamed)
    b.pointer_click('Array.from(document.querySelectorAll("[role=dialog] button")).find(e=>e.textContent==="Apply document settings")')
    b.wait('!document.querySelector("[role=dialog]")')
    selected(b, chosen['id'], renamed)
    b.wait('(' + row_expression + ').querySelector("strong").textContent===' + json.dumps(renamed))
    b.wait('(' + tab_expression + ').querySelector("span").textContent===' + json.dumps(renamed))
    assert b.evaluate('(' + row_expression + ').dataset.documentName') == renamed
    assert b.evaluate('(' + row_expression + ').dataset.documentLength') == old_length
    assert b.evaluate('(' + row_expression + ').querySelector("small").textContent') == old_metadata
    for expression in (row_expression, tab_expression):
        assert b.evaluate('(' + expression + ').title').startswith(renamed + '\n' + old_metadata + '\nDocument ID: ' + chosen['id'])
    deadline = time.monotonic() + 30
    while True:
        after_edit = persisted(b)
        edited_doc = next(doc for doc in after_edit[0]['documents'] if doc['id'] == chosen['id'])
        if edited_doc['head'] != previous_doc['head']:
            break
        assert time.monotonic() < deadline, 'Renamed revision did not persist'
        time.sleep(.1)
    assert [doc for doc in before_edit[0]['documents'] if doc['id'] != chosen['id']] == [doc for doc in after_edit[0]['documents'] if doc['id'] != chosen['id']]
    old_events = {event['id']: event for event in before_edit[0]['history']}
    new_events = [event for event in after_edit[0]['history'] if event['id'] not in old_events]
    assert len(new_events) == 1 and new_events[0]['document_id'] == chosen['id'], new_events
    assert new_events[0]['parents'] == [previous_doc['head_event']]
    assert new_events[0]['revision'] == edited_doc['head']
    assert all(event in after_edit[0]['history'] for event in old_events.values())
    assert all(after_edit[1].get(name) == value for name, value in before_edit[1].items()), 'Rename discarded historical source blobs'
    b.call('Page.reload')
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    selected(b, chosen['id'], renamed)
    assert_unchanged(after_edit, persisted(b))
    viewport(b, 1200)
    pagination(b, after_edit)
    assert not b.exceptions and not b.console_errors, (b.exceptions, b.console_errors)
    print('PASS: real annotated GenBank imports; compact searchable/sortable document explorer; row/tab identity switching; keyboard arrows/Home/End and wrapping; rail collapse preserves selection; 1600/1280/1200 desktop layout; native source/history/revision preservation and reload; same-ID reactive rename and exact one-event history; 101-record pagination/tab reveal and source preservation; no fresh runtime errors', flush=True)


if __name__ == '__main__':
    main()
