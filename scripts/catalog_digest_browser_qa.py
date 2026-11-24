#!/usr/bin/env python3
"""Desktop acceptance for pinned catalog geometry, documentary evidence and fresh products.

No scientific library is an oracle: motifs/cut offsets and sequence expectations
below are independent literal fixtures. Mobile review remains deferred.
"""
import io
import json
import re
import zipfile
from pathlib import Path

from browser_qa import Browser
from methods_browser_qa import fill
from workspace_browser_qa import slots
import desktop_mvp_browser_qa as mvp

DIRECTORY = Path('/tmp/collagene2-catalog-digest-qa')
ORIGIN = 'http://127.0.0.1:8080'
STARTER = ['EcoRI', 'BamHI', 'HindIII', 'EcoRV', 'SmaI', 'NotI', 'XbaI', 'BsaI', 'PstI', 'BsmBI', 'BbsI', 'SapI']
GEOMETRY = {
    'XhoI': ('CTCGAG', 1, 5), 'HpaII': ('CCGG', 1, 3), 'MspI': ('CCGG', 1, 3),
    'SmaI': ('CCCGGG', 3, 3), 'XmaI': ('CCCGGG', 1, 5),
    'Acc65I': ('GGTACC', 1, 5), 'KpnI': ('GGTACC', 5, 1), 'HinfI': ('GANTC', 1, 4),
}
MASK = {'A': set('A'), 'C': set('C'), 'G': set('G'), 'T': set('T'), 'N': set('ACGT')}


def set_enzymes(browser, selected):
    fill(browser, 'Find digest enzyme', '')
    labels = browser.evaluate('Array.from(document.querySelectorAll(".digest-enzyme-choice input")).map(e=>e.dataset.enzyme)')
    assert len(labels) == 40 and len(set(labels)) == 40
    for label in labels:
        target = 'document.querySelector(' + json.dumps('[aria-label="Digest enzyme ' + label + '"]') + ')'
        if browser.evaluate('(' + target + ').checked') != (label in selected):
            browser.pointer_click(target)
    actual = browser.evaluate('Array.from(document.querySelectorAll(".digest-enzyme-choice input:checked")).map(e=>e.dataset.enzyme)')
    assert set(actual) == set(selected), (actual, selected)


def digest(browser, selected):
    mvp.click(browser, 'Digest')
    set_enzymes(browser, selected)
    mvp.click(browser, 'Run virtual digest')
    browser.wait('!!document.querySelector(".digest-results")', timeout=60)
    return json.loads(mvp.download(browser, 'Download digest evidence', 'virtual-digest.json'))


def literal_sites(sequence, circular, enzyme):
    motif, top, bottom = GEOMETRY[enzyme]
    length = len(sequence)
    found = []
    # Every motif here is reverse-complement symmetric, including IUPAC N.
    for start in range(length if circular else max(0, length - len(motif) + 1)):
        observed = ''.join(sequence[(start + i) % length] for i in range(len(motif)))
        if all(base in MASK[code] for base, code in zip(observed, motif)):
            t, b = start + top, start + bottom
            found.append((start, t % length if circular else t, b % length if circular else b))
    return found


def assert_geometry(artifact, sequence, circular, selected):
    report = artifact['report']
    assert report['scan_complete']
    for enzyme in selected:
        observed = [(s['recognition_start'], s['top_cut'], s['bottom_cut']) for s in report['sites'] if s['enzyme'] == enzyme]
        assert observed == literal_sites(sequence, circular, enzyme), (enzyme, observed)
    selection = artifact['catalog_pin']
    assert selection['schema_version'] == 1
    assert re.fullmatch('[0-9a-f]{64}', selection['catalog_content_hash'])
    assert artifact['catalog_engine_version'] == 'catalog-digest-1'
    assert re.fullmatch('[0-9a-f]{64}', artifact['catalog_batch_revision'])
    assert [e['enzyme'] for e in selection['enzymes']] == report['enzymes']
    assert {e['enzyme']['name'] for e in selection['enzymes']} == set(selected)


def reveal_catalog(browser):
    if not browser.evaluate('document.querySelector(".digest-catalog").open'):
        browser.pointer_click('document.querySelector(".digest-catalog > summary")')


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for path in DIRECTORY.iterdir():
        if path.is_file():
            path.unlink()
    mvp.DIRECTORY = DIRECTORY
    source_sequence = 'AAAACTCGAG' + 'T' * 20 + 'AAAA'
    relations_sequence = 'AAAACCCGGGAAAAGGTACCAAAACCGGAAAA'
    wrap_sequence = 'TC' + 'A' * 17 + 'GAC'
    fixtures = [
        ('CATALOG_SOURCE', source_sequence, True, '     misc_feature    13..22\n                     /label="Catalog retained annotation"\n'),
        ('CATALOG_RELATIONS', relations_sequence, False, ''),
        ('CATALOG_WRAP', wrap_sequence, True, ''),
        ('CATALOG_UNCUT', 'A' * 47, True, ''),
    ]
    paths = {}
    for name, sequence, circular, features in fixtures:
        path = DIRECTORY / (name + '.gb')
        path.write_text(mvp.record(name, sequence, circular, features))
        paths[name] = path
    browser = Browser()
    browser.call('Browser.setDownloadBehavior', {'behavior': 'allow', 'downloadPath': str(DIRECTORY)})
    browser.call('Emulation.setDeviceMetricsOverride', {'width': 1600, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})
    browser.call('Storage.clearDataForOrigin', {'origin': ORIGIN, 'storageTypes': 'indexeddb'})
    # Clear only historical Runtime.enable replay, before navigating the fresh app.
    browser.exceptions.clear()
    browser.console_errors.clear()
    browser.call('Page.navigate', {'url': ORIGIN})
    browser.wait('document.body?.innerText.includes("Saved in this browser")')
    mvp.click(browser, 'Digest')
    labels = browser.evaluate('Array.from(document.querySelectorAll(".digest-enzyme-choice input")).map(e=>e.dataset.enzyme)')
    assert labels[:12] == STARTER and len(labels) == 40
    assert browser.evaluate('Array.from(document.querySelectorAll(".digest-enzyme-choice input:checked")).map(e=>e.dataset.enzyme)') == ['EcoRI']
    assert browser.evaluate('Array.from(document.querySelectorAll(".mvp-inspector-tabs > button")).map(e=>e.textContent).join(",")') == 'Features,Digest'
    fill(browser, 'Find digest enzyme', 'MspI')
    assert browser.evaluate('Array.from(document.querySelectorAll(".digest-enzyme-choice input")).map(e=>e.dataset.enzyme)') == ['MspI']
    fill(browser, 'Find digest enzyme', 'CTCGAG')
    assert browser.evaluate('Array.from(document.querySelectorAll(".digest-enzyme-choice input")).map(e=>e.dataset.enzyme)') == ['XhoI']
    fill(browser, 'Find digest enzyme', '')
    for name, *_ in fixtures:
        mvp.upload(browser, paths[name])
        browser.wait('document.querySelector("h1")?.textContent===' + json.dumps(name))

    # A newly available geometry drives the same fresh product and history path.
    mvp.document(browser, 'CATALOG_SOURCE')
    evidence = digest(browser, ['XhoI'])
    assert_geometry(evidence, source_sequence, True, ['XhoI'])
    selection = evidence['catalog_pin']
    assert [e['id'] for e in selection['enzymes']] == ['xhoi']
    assert selection['catalog_edition'] == browser.evaluate('document.querySelector(".digest-catalog").dataset.catalogEdition')
    assert selection['catalog_content_hash'] == browser.evaluate('document.querySelector(".digest-catalog").dataset.catalogHash')
    report = evidence['report']
    assert report['products_determined'] and report['ref_top_lengths'] == [len(source_sequence)]
    assert report['fragments'][0]['top_start'] == 5
    assert report['fragments'][0]['duplex']['ends'][0]['overhang'] == 'TCGA'
    reveal_catalog(browser)
    assert 'Unknown:' in browser.evaluate('document.querySelector(".digest-catalog-identity").textContent')
    exported = mvp.download(browser, 'Export fragment GenBank', 'digest-fragment-1.gb').decode()
    assert mvp.gb_sequence(exported) == source_sequence[5:] + source_sequence[:5]
    assert 'Catalog retained annotation' in exported and '8..17' in exported
    assert mvp.digest_comment(exported, 'CG-CATALOG-DIGEST-v1') == selection
    assert mvp.digest_comment(exported, 'CG-DIGEST-RECIPE-v1')['inspection_revision'] == report['batch_revision']
    assert 'unverified' in exported.lower()

    # Catalog products use the existing atomic transaction, including real failure/abort.
    browser.wait('document.body?.innerText.includes("Saved in this browser")')
    previous_slots = slots(browser)
    count = browser.evaluate('document.querySelectorAll(".document-entry").length')
    browser.evaluate('window.__put=IDBObjectStore.prototype.put;IDBObjectStore.prototype.put=function(value,key){if(this.transaction.db.name==="collagene"&&this.name==="workspace"&&key==="current")throw new DOMException("Catalog quota fixture","QuotaExceededError");return window.__put.call(this,value,key)}')
    mvp.click(browser, 'Open fragment')
    browser.wait('document.querySelector(".digest-open-panel")?.textContent.includes("fragment was not opened")', timeout=60)
    browser.evaluate('IDBObjectStore.prototype.put=window.__put')
    assert slots(browser) == previous_slots
    assert browser.evaluate('document.querySelectorAll(".document-entry").length') == count
    browser.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()')
    browser.wait('!document.querySelector("[role=dialog]")')
    browser.evaluate('window.__put=IDBObjectStore.prototype.put;window.__cancelledPending=false;IDBObjectStore.prototype.put=function(value,key){const r=window.__put.call(this,value,key);if(this.transaction.db.name==="collagene"&&this.name==="workspace"&&key==="current"){window.__cancelledPending=true;queueMicrotask(()=>document.querySelector("button[aria-label=\\"Close dialog\\"]")?.click())}return r}')
    mvp.click(browser, 'Open fragment')
    browser.wait('window.__cancelledPending&&!document.querySelector("[role=dialog]")&&!document.querySelector(".application").inert', timeout=60)
    browser.evaluate('IDBObjectStore.prototype.put=window.__put')
    assert slots(browser) == previous_slots
    assert browser.evaluate('document.querySelectorAll(".document-entry").length') == count
    mvp.click(browser, 'Open fragment')
    browser.wait('!document.querySelector("[role=dialog]")&&document.querySelector("h1")?.textContent==="CATALOG_SOURCE_digest_fragment_1"', timeout=60)
    workspace = mvp.archive(browser)
    with zipfile.ZipFile(io.BytesIO(workspace)) as archive:
        manifest = json.loads(archive.read('manifest.json'))
        product = next(d for d in manifest['documents'] if d['id'] == manifest['active_document'])
        event = next(e for e in manifest['history'] if e['id'] == product['head_event'])
        assert event['parents'] == [evidence['source']['binding']['event_id']]
    native = mvp.download(browser, 'Save .cg2', 'CATALOG_SOURCE_digest_fragment_1.cg2')
    assert b'restriction_catalog:selection' in native
    mvp.click(browser, 'Digest')
    assert not browser.evaluate('!!document.querySelector(".digest-results")')
    retained = json.loads(mvp.download(browser, 'Download retained catalog selection', 'digest-retained-catalog.json'))
    assert retained == selection
    assert 'unverified' in browser.evaluate('document.querySelector(".digest-provenance").textContent').lower()
    browser.call('Page.reload')
    browser.wait('document.body?.innerText.includes("Saved in this browser")')
    mvp.click(browser, 'Digest')
    assert not browser.evaluate('!!document.querySelector(".digest-results")')
    assert json.loads(mvp.download(browser, 'Download retained catalog selection', 'digest-retained-catalog.json')) == selection
    # Text import retains explicitly unverified descriptive catalog claims too.
    mvp.upload(browser, DIRECTORY / 'digest-fragment-1.gb')
    browser.wait('document.querySelectorAll(".document-entry").length>' + str(count + 1))
    mvp.click(browser, 'Digest')
    assert not browser.evaluate('!!document.querySelector(".digest-results")')
    assert json.loads(mvp.download(browser, 'Download retained catalog selection', 'digest-retained-catalog.json')) == selection

    # Isoschizomers stay distinct, and neoschizomers retain different paired cuts.
    mvp.document(browser, 'CATALOG_RELATIONS')
    for selected in [['HpaII', 'MspI'], ['SmaI'], ['XmaI'], ['Acc65I'], ['KpnI']]:
        result = digest(browser, selected)
        assert_geometry(result, relations_sequence, False, selected)
        if selected == ['HpaII', 'MspI']:
            assert len({e['id'] for e in result['catalog_pin']['enzymes']}) == 2
            by_enzyme = {name: [(s['recognition_start'], s['top_cut'], s['bottom_cut']) for s in result['report']['sites'] if s['enzyme'] == name] for name in selected}
            assert by_enzyme['HpaII'] == by_enzyme['MspI'] and by_enzyme['HpaII']
    # Documentary facts are inspectable, formulation-specific, and explicitly noncausal.
    digest(browser, ['EcoRI', 'SmaI'])
    reveal_catalog(browser)
    for ident in ['neb-r0101', 'neb-r3101', 'neb-r0141']:
        target = 'document.querySelector(' + json.dumps('[data-formulation-id="' + ident + '"] > summary') + ')'
        browser.pointer_click(target)
    text = browser.evaluate('document.querySelector(".digest-catalog").innerText')
    assert '50% relative activity' in text and '100% relative activity' in text
    assert 'Context:' in text and 'star activity' in text
    assert '25 C' in text and 'longer than one hour' in text
    assert 'do not establish cleavage' in text
    assert browser.evaluate('Array.from(document.querySelectorAll(".digest-catalog-citations a")).every(a=>a.href.startsWith("https://"))')
    assert browser.evaluate('!!document.querySelector(".digest-catalog-citations a[href*=which-incubation-temperature-should-i-use]")')
    browser.evaluate('document.querySelector("[data-formulation-id=neb-r0141] .digest-catalog-fact").scrollIntoView({block:"start"})')
    browser.screenshot('/tmp/collagene2-catalog-digest-evidence.png')

    # IUPAC recognition spans the circular origin with its real five-base footprint.
    mvp.document(browser, 'CATALOG_WRAP')
    wrapped = digest(browser, ['HinfI'])
    assert_geometry(wrapped, wrap_sequence, True, ['HinfI'])
    assert wrapped['report']['products_determined']
    assert {end['overhang'] for end in wrapped['report']['fragments'][0]['duplex']['ends']} == {'ACT', 'AGT'}
    site = wrapped['report']['sites'][0]
    assert site['recognition_length'] == 5 and sum(s['end'] - s['start'] for s in site['spans']) == 5 and len(site['spans']) == 2
    mvp.click_expression(browser, 'document.querySelector(".digest-site button")')
    assert browser.evaluate('document.querySelector(".statusbar").textContent.includes("Click a base")')
    # Every catalog identity is supported by the same bounded inspection path.
    mvp.document(browser, 'CATALOG_UNCUT')
    all_names = browser.evaluate('Array.from(document.querySelectorAll(".digest-enzyme-choice input")).map(e=>e.dataset.enzyme)')
    all_result = digest(browser, all_names)
    assert len(all_result['catalog_pin']['enzymes']) == 40 and len(all_result['report']['enzymes']) == 40
    assert all_result['report']['products_determined'] and all_result['report']['ref_top_lengths'] == [47]
    assert not all_result['report']['sites']
    # Resource limits and source/options cancellation cannot restore stale results.
    fill(browser, 'Digest comparison limit', '1')
    mvp.click(browser, 'Run virtual digest')
    browser.wait('!!document.querySelector(".digest-results")')
    capped = json.loads(mvp.download(browser, 'Download digest evidence', 'virtual-digest.json'))
    assert not capped['report']['scan_complete'] and not capped['report']['products_determined']
    fill(browser, 'Digest comparison limit', '50000000')
    mvp.hold_jobs(browser)
    mvp.click(browser, 'Run virtual digest')
    browser.wait('window.__heldPosts.length>0')
    mvp.click(browser, 'Cancel digest')
    mvp.release_jobs(browser)
    assert not browser.evaluate('!!document.querySelector(".digest-results")')
    mvp.hold_jobs(browser)
    mvp.click(browser, 'Run virtual digest')
    browser.wait('window.__heldPosts.length>0')
    set_enzymes(browser, ['XhoI'])
    mvp.release_jobs(browser)
    assert not browser.evaluate('!!document.querySelector(".digest-results")')
    mvp.document(browser, 'CATALOG_SOURCE')
    original = digest(browser, ['XhoI'])
    mvp.upload(browser, paths['CATALOG_SOURCE'])
    browser.wait('Array.from(document.querySelectorAll(".document-entry strong")).filter(e=>e.textContent==="CATALOG_SOURCE").length===2')
    mvp.document(browser, 'CATALOG_SOURCE', 0)
    mvp.hold_jobs(browser)
    mvp.click(browser, 'Run virtual digest')
    browser.wait('window.__heldPosts.length>0')
    browser.evaluate('window.__holdJobs=false')
    mvp.document(browser, 'CATALOG_SOURCE', 1)
    mvp.release_jobs(browser)
    assert not browser.evaluate('!!document.querySelector(".digest-results")')
    duplicate = digest(browser, ['XhoI'])
    assert duplicate['source']['binding']['revision'] == original['source']['binding']['revision']
    assert duplicate['source']['binding']['document_id'] != original['source']['binding']['document_id']
    assert duplicate['catalog_batch_revision'] != original['catalog_batch_revision']
    for width in [1600, 1280, 1200]:
        browser.call('Emulation.setDeviceMetricsOverride', {'width': width, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})
        browser.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')
        metrics = browser.evaluate('({width:innerWidth,document:document.documentElement.scrollWidth,panel:document.querySelector(".digest-panel").scrollWidth,client:document.querySelector(".digest-panel").clientWidth})')
        assert metrics['document'] <= width and metrics['panel'] <= metrics['client'] + 1, metrics
    browser.screenshot('/tmp/collagene2-catalog-digest-desktop.png')
    assert not browser.exceptions and not browser.console_errors, (browser.exceptions, browser.console_errors)
    print('PASS: catalog40/original12 controls, search, exact pinned geometry and ordered identities, independent neo/isoschizomer/IUPAC circular oracles, formulation contexts/qualifiers/citations/unknowns, fresh annotated GenBank and catalog provenance, atomic native install/rollback/abort with exact source parent, unverified native/text/reload evidence, resource caps/cancellation/identical-content source guards, desktop layout')


if __name__ == '__main__':
    main()
