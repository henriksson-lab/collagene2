#!/usr/bin/env python3
"""Development-only annotation, translation and pairwise UI regression checks."""
import json
from pathlib import Path

from browser_qa import QA_ORIGIN, Browser, mobile_enabled
from format_browser_qa import wait_file


def main():
    b = Browser()
    directory = Path('/tmp/collagene2-biology-qa')
    directory.mkdir(exist_ok=True)
    for path in directory.iterdir():
        if path.is_file():
            path.unlink()
    b.call('Browser.setDownloadBehavior', {'behavior': 'allow', 'downloadPath': str(directory)})
    b.call('Emulation.setDeviceMetricsOverride', {'width': 1600, 'height': 1000, 'deviceScaleFactor': 1, 'mobile': False})
    b.call('Storage.clearDataForOrigin', {'origin': QA_ORIGIN, 'storageTypes': 'indexeddb'})
    b.call('Page.navigate', {'url': QA_ORIGIN})
    b.wait("document.querySelector('h1')?.textContent==='pCG_demo'")
    b.wait("document.body?.innerText.includes('Saved in this browser')")

    def field(label, value, event='input'):
        selector = '[aria-label=' + json.dumps(label) + ']'
        b.evaluate('(()=>{const e=document.querySelector(' + json.dumps(selector) + ');if(!e)throw Error("Missing field");e.value=' + json.dumps(value) + ';e.dispatchEvent(new Event(' + json.dumps(event) + ',{bubbles:true}));})()')

    def paste(text, name):
        b.click('Paste sequence')
        field('Sequence input', text)
        b.click('Open sequence')
        b.wait('document.querySelector("h1")?.textContent===' + json.dumps(name))

    def feature(label):
        b.evaluate('Array.from(document.querySelectorAll(".feature-entry")).find(b=>b.querySelector("strong").textContent===' + json.dumps(label) + ').click()')

    def undo():
        b.evaluate('document.querySelector("button[title=Undo]").click()')

    paste('>biology_fixture\nATGAAATAACCC\n', 'biology_fixture')
    b.click('Biology')
    b.click('Translate frame')
    b.wait("document.querySelector('.protein-preview')?.textContent==='MK*P'")
    b.click('Download protein FASTA')
    assert b'MK*P' in wait_file(directory / 'translation.fasta')
    field('Translation frame', '1', 'change')
    b.wait("!document.querySelector('.protein-preview')")
    field('Translation frame', '0', 'change')
    field('ORF minimum amino acids', '1')
    b.evaluate("""(()=>{const e=document.querySelector('[aria-label="ORF both strands"]');e.checked=false;e.dispatchEvent(new Event("change",{bubbles:true}));})()""")
    b.click('Find ORFs')
    b.wait("document.querySelector('.orf-result h3')?.textContent==='1 ORF candidates'")
    b.click('Annotate ORF 1')
    b.click('Features')
    b.wait("document.querySelector('.feature-list').textContent.includes('Predicted ORF 1')")
    feature('Predicted ORF 1')
    assert b.evaluate("document.querySelector('.annotation-sequence').textContent==='ATGAAATAA'")
    assert b.evaluate("document.querySelector('.feature-entry').textContent.includes('review required')")
    b.click('Biology')
    b.click('Translate selected CDS')
    b.wait("document.querySelector('.biology-status')?.textContent.includes('requires review')")
    b.click('Features')
    b.click('Edit annotation')
    b.evaluate("""(()=>{const e=document.querySelector('.feature-editor input[type=checkbox]');e.checked=false;e.dispatchEvent(new Event('change',{bubbles:true}));})()""")
    b.click('Save annotation')
    feature('Predicted ORF 1')
    b.click('Biology')
    b.click('Translate selected CDS')
    b.wait("document.querySelector('.protein-preview')?.textContent==='MK'")

    # Independent compound location; exact extraction is reverse(complement(join)).
    b.click('Features')
    b.click('New annotation')
    field('Annotation name', 'compound_reverse')
    field('Annotation type', 'misc_feature')
    field('Annotation location', 'complement(join(10..12,1..3))')
    b.click('Add qualifier')
    b.evaluate("""(()=>{const e=document.querySelector('[aria-label="Qualifier value"]');e.value="first evidence";e.dispatchEvent(new Event("input",{bubbles:true}));})()""")
    b.click('Add qualifier')
    b.evaluate("""(()=>{const e=document.querySelectorAll('[aria-label="Qualifier value"]')[1];e.value="second evidence";e.dispatchEvent(new Event("input",{bubbles:true}));})()""")
    b.click('Save annotation')
    feature('compound_reverse')
    b.wait("document.querySelector('.annotation-sequence')?.textContent==='CATGGG'")
    assert b.evaluate("Array.from(document.querySelectorAll('button')).find(b=>b.textContent.trim()==='Edit selection').disabled")
    assert b.evaluate("document.querySelector('.sr-only[aria-live]').textContent.includes('CATGGG')")
    b.click('Edit annotation')
    assert b.evaluate("Array.from(document.querySelectorAll('[aria-label=\"Qualifier value\"]')).map(e=>e.value).join('|')==='first evidence|second evidence'")
    field('Annotation name', 'renamed_compound')
    b.click('Save annotation')
    feature('renamed_compound')
    b.click('Delete annotation')
    b.wait("!document.querySelector('.feature-list').textContent.includes('renamed_compound')")
    undo()
    b.wait("document.querySelector('.feature-list').textContent.includes('renamed_compound')")
    b.click('Document settings')
    field('Document name', 'renamed_document')
    field('Reference topology', 'circular', 'change')
    b.click('Apply document settings')
    b.wait("document.querySelector('h1')?.textContent==='renamed_document'")
    b.click('Document settings')
    field('New origin', '4')
    b.click('Rotate origin')
    b.click('FASTA')
    assert b'AAATAACCCATG' in wait_file(directory / 'renamed_document.fasta')
    (directory / 'renamed_document.fasta').unlink()
    undo()
    b.click('Document settings')
    b.click('Reverse complement document')
    b.click('FASTA')
    assert b'GGGTTATTTCAT' in wait_file(directory / 'renamed_document.fasta')
    undo()

    # Separate reference/query documents test semiglobal coordinates and reactive clearing.
    paste('>pair_reference\nTTTACGTGGG\n>pair_query\nACGT\n', 'pair_reference')
    b.click('Compare')
    # Query indices follow the current open-document order.
    query_index = b.evaluate("Array.from(document.querySelector('[aria-label=\"Query document\"]').options).find(o=>o.textContent==='pair_query').value")
    field('Query document', query_index, 'change')
    field('Alignment mode', 'query', 'change')
    b.click('Compare sequences')
    b.wait("document.querySelector('.compare-result h3')?.textContent==='Score 8'")
    assert b.evaluate("document.querySelector('.compare-result').textContent.includes('Reference [3, 7); query [0, 4).')")
    assert b.evaluate("document.querySelector('.compare-result').textContent.includes('reference 40.0% · query 100.0%')")
    b.click('Download comparison JSON')
    report = json.loads(wait_file(directory / 'pairwise-comparison.json'))
    assert report['alignment']['cigar'] == '4='
    assert report['reference_revision'] != report['query_revision']
    field('Mismatch', '-5')
    b.wait("!document.querySelector('.compare-result')")
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride', {'width': 390, 'height': 844, 'deviceScaleFactor': 1, 'mobile': True})
        b.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        assert b.evaluate('document.documentElement.scrollWidth<=window.innerWidth'), 'Mobile overflow'
    assert not b.exceptions, b.exceptions
    assert not b.console_errors, b.console_errors
    if mobile_enabled():
        b.screenshot('/tmp/collagene2-biology-mobile.png')
    print('PASS: annotation create/edit/compound/qualifiers/delete/undo; topology/origin/reverse; translation/CDS/ORF; pairwise worker/export/stale options; mobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__ == '__main__':
    main()
