#!/usr/bin/env python3
"""Development-only GFF/BED mapping, feature library, review and evidence QA."""
import json
import time
from pathlib import Path
from browser_qa import QA_ORIGIN, Browser, mobile_enabled


def field(b,label,value,event='input'):
    selector=json.dumps('[aria-label='+json.dumps(label)+']')
    b.evaluate(f'(()=>{{const e=document.querySelector({selector});if(!e)throw Error("Missing field {label}");e.value={json.dumps(value)};e.dispatchEvent(new Event({json.dumps(event)},{{bubbles:true}}));}})()')


def paste(b,text):
    b.click('Paste sequence');field(b,'Sequence input',text);b.click('Open sequence')


def file_bytes(path):
    until=time.monotonic()+10
    while not path.exists() and time.monotonic()<until:time.sleep(.1)
    assert path.exists(),path
    return path.read_bytes()


def choose(b,name):
    b.evaluate('Array.from(document.querySelectorAll(".document-entry")).find(e=>e.querySelector("strong").textContent==='+json.dumps(name)+').click()')
    b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name))


def main():
    b=Browser();directory=Path('/tmp/collagene2-annotation-library-qa');directory.mkdir(exist_ok=True)
    for path in directory.iterdir():
        if path.is_file():path.unlink()
    b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(directory)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':QA_ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':QA_ORIGIN})
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    paste(b,'>annotation_source\nTTTTACGTTGGG\n')
    b.wait("document.querySelector('h1')?.textContent==='annotation_source'")
    b.click('Annotation files')
    gff='##gff-version 3\n##sequence-region refA 1 12\nrefA\tfixture\tpromoter\t5\t9\t.\t+\t.\tID=motif;Name=My%20motif;Note=Retained%20source\nrefB\tfixture\tgene\t1\t3\t.\t+\t.\tID=other\n'
    field(b,'Annotation reference ID','refA');field(b,'Annotation source text',gff)
    b.click('Review annotation import')
    b.wait("document.querySelector('.annotation-file-report')?.textContent.includes('1 imported annotation proposals')")
    assert b.evaluate("document.querySelector('.annotation-file-report').textContent.includes('refB')")
    assert b.evaluate("document.querySelector('.annotation-file-proposal').textContent.includes('5..9')")
    b.click('Download original annotation source');assert file_bytes(directory/'original-annotations.txt')==gff.encode()
    b.click('Select imported annotations on page');b.click('Apply selected imported annotations')
    b.wait("!document.querySelector('.annotation-file-report')")
    b.click('Features');b.wait("document.querySelector('.feature-list')?.textContent.includes('My motif')")
    b.evaluate("document.querySelector('.feature-entry').click()")
    b.wait("document.querySelector('.annotation-sequence')?.textContent==='ACGTT'")
    b.click('Use selected annotation in library')
    field(b,'Feature library name','My motifs');field(b,'Feature library version','v1')
    b.click('Create feature library');b.click('Add selected annotation to library')
    b.wait("document.querySelector('.template-list')?.textContent.includes('My motif · 5 bases')")
    b.click('Download feature library')
    lib=json.loads(file_bytes(directory/'feature-library.cg2features'))
    assert lib['entries'][0]['sequence']['bases']=='ACGTT',lib['entries'][0]['sequence']
    paste(b,'>library_target\nGGGACGTTCCCAACGTAAA\n')
    b.wait("document.querySelector('h1')?.textContent==='library_target'")
    b.click('Feature libraries')
    # The selected immutable library revision may reset when the panel is reopened.
    value=b.evaluate("Array.from(document.querySelector('[aria-label=\"Feature library revision\"]').options).find(e=>e.textContent.includes('1 entries')).value")
    field(b,'Feature library revision',value,'change')
    b.click('Search feature library');b.wait("document.querySelector('.feature-library-report')?.textContent.includes('2 annotation proposals')")
    assert b.evaluate("Array.from(document.querySelectorAll('.library-hit')).every(e=>e.textContent.includes('100.00–100.00%'))")
    b.evaluate("document.querySelector('.feature-library-report').scrollIntoView({block:'start'})")
    b.screenshot('/tmp/collagene2-library-proposals.png')
    b.evaluate("document.querySelectorAll('.library-hit input[type=checkbox]').forEach(e=>{e.checked=true;e.dispatchEvent(new Event('change',{bubbles:true}));})")
    b.click('Apply selected library annotations');b.wait("!document.querySelector('.feature-library-report')")
    b.click('Features');b.wait("document.querySelectorAll('.feature-entry').length===2")
    assert b.evaluate("Array.from(document.querySelectorAll('.feature-entry')).every(e=>e.textContent.includes('review required'))")
    # Independent source strings verify reverse biological orientation.
    b.evaluate("document.querySelectorAll('.feature-entry')[1].click()")
    b.wait("document.querySelector('.annotation-sequence')?.textContent==='ACGTT'")
    b.evaluate("document.querySelector('button[title=Undo]').click()")
    b.wait("document.querySelectorAll('.feature-entry').length===0")
    b.evaluate("document.querySelector('button[title=Redo]').click()")
    b.wait("document.querySelectorAll('.feature-entry').length===2")
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    b.call('Page.reload',{'ignoreCache':True});b.wait("document.body?.innerText.includes('Saved in this browser')")
    b.wait("document.querySelector('h1')?.textContent==='library_target'")
    b.click('Feature libraries')
    b.wait("document.querySelector('[aria-label=\"Feature library revision\"]')?.textContent.includes('1 entries')")
    value=b.evaluate("Array.from(document.querySelector('[aria-label=\"Feature library revision\"]').options).find(e=>e.textContent.includes('1 entries')).value")
    field(b,'Feature library revision',value,'change');b.click('Search feature library')
    b.wait("document.querySelectorAll('.library-hit').length===2")
    assert b.evaluate("Array.from(document.querySelectorAll('.library-hit input')).every(e=>e.disabled)")
    # BED12 is a separate explicit mapping; preserve compound reverse biological order.
    b.click('Annotation files');field(b,'Annotation file format','bed','change');field(b,'Annotation reference ID','chrTarget')
    bed='track name="Fixture"\nchrTarget\t0\t19\tBlocked feature\t0\t-\t0\t0\t0\t2\t3,3\t0,16\n'
    field(b,'Annotation source text',bed);b.click('Review annotation import')
    b.wait("document.querySelector('.annotation-file-report')?.textContent.includes('1 imported annotation proposals')")
    assert b.evaluate("document.querySelector('.annotation-file-proposal').textContent.includes('complement(join(1..3,17..19))')")
    b.click('Select imported annotations on page');b.click('Apply selected imported annotations')
    b.wait("!document.querySelector('.annotation-file-report')")
    b.click('Features');b.wait("document.querySelectorAll('.feature-entry').length===3")
    b.evaluate("document.querySelectorAll('.feature-entry')[2].click()")
    b.wait("document.querySelector('.annotation-sequence')?.textContent==='TTTCCC'")
    b.click('Annotation files');field(b,'Annotation reference ID','chrTarget');field(b,'Annotation file format','bed','change');b.click('Export current annotations')
    assert b'chrTarget\t' in file_bytes(directory/'annotations.bed')
    choose(b,'annotation_source');b.click('Annotation files')
    b.wait("!!document.querySelector('.retained-annotation-sources')")
    b.evaluate("document.querySelector('.retained-annotation-sources button').click()")
    assert file_bytes(directory/'retained-annotations.txt')==gff.encode()
    invalid=directory/'invalid-encoding.gff3';invalid.write_bytes(b'##gff-version 3\nref\tbad\tgene\t1\t2\t.\t+\t.\tName=\xff\n')
    before=b.evaluate("document.querySelector('[aria-label=\"Annotation source text\"]').value")
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId']
    node=b.call('DOM.querySelector',{'nodeId':root,'selector':'.annotation-files-panel input[type=file]'})['nodeId']
    b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(invalid)]})
    b.wait("document.querySelector('.annotation-file-status')?.textContent.includes('not valid UTF-8')")
    assert b.evaluate("document.querySelector('[aria-label=\"Annotation source text\"]').value")==before
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    b.screenshot('/tmp/collagene2-annotation-files.png')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        b.wait('document.documentElement.scrollWidth<=innerWidth')
    print('PASS: GFF3/BED coordinates/reference mapping/source evidence; reusable library orientation/review/apply/undo/reload/duplicate guard; worker annotation exports')
    b.socket.close()


if __name__=='__main__':main()
