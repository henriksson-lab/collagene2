"""Existing alignment exchange, column mapping and transactional persistence in WASM.
This creates only synthetic fixtures and clears the isolated test origin.
"""
import io
import json
import time
import zipfile
from pathlib import Path
from browser_qa import QA_ORIGIN, Browser, mobile_enabled
from workspace_browser_qa import wait_file

ROOT=Path('/tmp/collagene2-alignment-qa')
STOCKHOLM='''# STOCKHOLM 1.0
#=GF ID literal_fixture
#=GS row_a DE original row evidence
row_a AC.GT
row_b A-AGT
gap_only -----
#=GR row_a PP 12.45
#=GC RF xx.xx
//
'''


def fill(b,selector,value,event='input'):
    b.evaluate('(()=>{const e=document.querySelector('+json.dumps(selector)+');e.value='+json.dumps(value)+';e.dispatchEvent(new Event('+json.dumps(event)+',{bubbles:true}));})()')


def open_file(b,path,selector='input[type=file]'):
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId']
    node=b.call('DOM.querySelector',{'nodeId':root,'selector':selector})['nodeId']
    assert node,(selector,b.evaluate('document.body.innerText'))
    b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})
    if selector=='input[type=file]':b.wait('document.querySelector("#dialog-title")?.textContent==="Alignments"')


def close(b):
    b.evaluate('document.querySelector('+json.dumps('button[aria-label="Close dialog"]')+').click()')
    b.wait('!document.querySelector("[role=dialog]")')


def key(b,name):
    b.evaluate('document.querySelector(".alignment-canvas").dispatchEvent(new KeyboardEvent("keydown",{key:'+json.dumps(name)+',bubbles:true}))')


def selected(b,name,checked=True):
    selector='input[aria-label='+json.dumps('Select alignment row '+name)+']'
    b.evaluate('(()=>{const e=document.querySelector('+json.dumps(selector)+');if(e.checked!=='+json.dumps(checked)+')e.click();})()')


def collection_count(b):
    return b.evaluate('document.querySelectorAll('+json.dumps('select[aria-label="Retained alignment"] option')+').length')


def imported(b,count=1):
    b.wait('document.querySelector(".alignment-status")?.textContent.includes("Alignment retained")')
    b.wait('document.querySelectorAll('+json.dumps('select[aria-label="Retained alignment"] option')+').length==='+str(count))
    b.wait('Number(document.querySelector(".alignment-canvas")?.dataset.paintedCells)>0')


def main():
    ROOT.mkdir(exist_ok=True)
    for p in ROOT.iterdir():
        if p.is_file():p.unlink()
    b=Browser();b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(ROOT)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':QA_ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':QA_ORIGIN})
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    source=ROOT/'literal.sto';source.write_text(STOCKHOLM)
    open_file(b,source)
    b.wait('document.querySelector("#dialog-title")?.textContent==="Alignments"')
    b.click('Import alignment')
    b.wait('document.querySelector(".alignment-status")?.textContent.includes("Choose the supplied")')
    assert not b.evaluate('!!document.querySelector(".alignment-canvas")')
    fill(b,'select[aria-label="Alignment supplied alphabet"]','Dna','change')
    b.click('Import alignment');imported(b)
    # Literal gap mappings are coordinates in each ungapped row, never references.
    key(b,'ArrowRight');key(b,'ArrowRight')
    b.wait('document.querySelector(".alignment-readout")?.textContent.includes("Row row_a · alignment column 3 · - · gap; no residue position")')
    assert 'GR PP=.' in b.evaluate('document.querySelector(".alignment-readout").textContent')
    assert 'GC RF=.' in b.evaluate('document.querySelector(".alignment-readout").textContent')
    key(b,'ArrowDown')
    b.wait('document.querySelector(".alignment-readout")?.textContent.includes("Row row_b · alignment column 3 · A · ungapped row position 2")')
    assert 'Counts are not alignment confidence' in b.evaluate('document.querySelector(".alignment-readout").textContent')
    b.click('Original alignment source')
    assert wait_file(ROOT/'original-alignment.sto')==source.read_bytes()
    b.click('Save alignment file');native=wait_file(ROOT/'alignment.cg2aln');doc=json.loads(native)
    assert doc['rows'][0]['aligned']=='AC-GT' and doc['rows'][1]['aligned']=='A-AGT'
    assert doc['source']['alphabet']=='Dna' and len(doc['source']['blake3'])==64
    # Export loss is reviewed before download and invalidated on options change.
    b.click('Prepare alignment export')
    b.wait('document.querySelector(".alignment-export-review")?.textContent.includes("FASTA omits all Stockholm")')
    b.click('Download alignment export')
    fasta=wait_file(ROOT/'alignment.fasta').decode();assert '>row_a\nAC-GT\n' in fasta and '>gap_only\n-----\n' in fasta
    fill(b,'select[aria-label="Alignment export format"]','Clustal','change')
    b.wait('!document.querySelector(".alignment-export-review")')
    # All-gap extraction fails atomically; valid rows open as separate linear docs.
    b.evaluate('document.querySelector(".alignment-rows").open=true')
    selected(b,'gap_only');selected(b,'row_a');before=b.evaluate('document.querySelectorAll(".document-entry").length')
    b.click('Open selected ungapped rows')
    b.wait('document.querySelector(".alignment-status")?.textContent.includes("at least one")')
    assert b.evaluate('document.querySelectorAll(".document-entry").length')==before
    b.click('Clear row selection');selected(b,'row_a');selected(b,'row_b')
    b.click('Open selected ungapped rows')
    b.wait('document.querySelector(".alignment-status")?.textContent.includes("Opened 2")')
    assert b.evaluate('document.querySelectorAll(".document-entry").length')==before+2
    b.screenshot('/tmp/collagene2-alignment.png')
    close(b);b.wait('document.querySelector("h1")?.textContent==="row_a"')
    b.click('FASTA');out=wait_file(ROOT/'row_a.fasta').decode();assert ''.join(out.splitlines()[1:])=='ACGT'
    assert b.evaluate('document.querySelector(".document-meta").textContent.includes("Linear")')
    b.wait('document.body.innerText.includes("Saved in this browser")')
    # Saved collection keeps source bytes through reload and workspace exchange.
    b.call('Page.reload');b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.click('Alignments');b.wait('!!document.querySelector(".alignment-canvas")')
    assert collection_count(b)==1
    close(b);open_file(b,ROOT/'alignment.cg2aln');b.wait('document.querySelector("#dialog-title")?.textContent==="Alignments"')
    b.click('Import alignment');imported(b)
    # Invalid saved source identity does not add or replace an existing alignment.
    close(b);bad=json.loads(native);bad['source']['blake3']='0'*64;bad_path=ROOT/'tampered.cg2aln';bad_path.write_text(json.dumps(bad))
    open_file(b,bad_path);b.click('Import alignment');b.wait('document.querySelector(".alignment-status")?.textContent.includes("source hash")')
    assert collection_count(b)==1
    close(b)
    # Native archives retain the independently validated typed alignment collection.
    b.click('Workspace');b.click('Save workspace file');archive=wait_file(ROOT/'collagene-workspace.cg2workspace')
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        manifest=json.loads(z.read('manifest.json'));assert any(a['kind']=='alignment-document-v1' for a in manifest['attachments'].values())
    # Original archive reopens without losing its independent collection evidence.
    open_file(b,ROOT/'collagene-workspace.cg2workspace','[role=dialog] input[type=file][accept*=cg2workspace]')
    b.wait('!document.querySelector("[role=dialog]")');b.wait('document.body.innerText.includes("Saved in this browser")')
    b.click('Alignments');b.wait('!!document.querySelector(".alignment-canvas")')
    # Pending jobs are abandoned when import options change before worker readiness.
    b.evaluate('document.querySelector(".alignment-import").open=true')
    if b.evaluate('Array.from(document.querySelectorAll("button")).some(b=>b.textContent.trim()==="Use pasted alignment instead")'):b.click('Use pasted alignment instead')
    fill(b,'textarea[aria-label="Aligned sequence input"]','>abandoned\nACGT\n')
    fill(b,'select[aria-label="Alignment source format"]','AlignedFasta','change')
    fill(b,'select[aria-label="Alignment supplied alphabet"]','Dna','change')
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(b=>b.textContent.trim()==="Import alignment").click();const e=document.querySelector('+json.dumps('input[aria-label="Alignment name"]')+');e.value="changed_during_import";e.dispatchEvent(new Event("input",{bubbles:true}));})()')
    time.sleep(.5)
    assert collection_count(b)==1
    close(b)
    # CLUSTAL file routing and imported consensus are separately preserved.
    clustal=ROOT/'literal.aln';clustal.write_text('CLUSTAL W synthetic\n\na AC-GT 4\nb A-AGT 4\n  *  **\n')
    open_file(b,clustal);fill(b,'select[aria-label="Alignment supplied alphabet"]','Dna','change');b.click('Import alignment');imported(b,2)
    b.evaluate('document.querySelector(".alignment-metadata").open=true')
    assert 'CLUSTAL consensus' in b.evaluate('document.querySelector(".alignment-metadata").textContent')
    close(b)
    # Large logical width is reachable with capped native CSS dimensions.
    wide=ROOT/'wide.afa';wide.write_text('>long\n'+'\n'.join(('A'*8000 for _ in range(500)))+'\n')
    open_file(b,wide);fill(b,'select[aria-label="Alignment supplied alphabet"]','Dna','change');b.click('Import alignment');imported(b,3)
    key(b,'End');b.wait('document.querySelector(".alignment-readout").textContent.includes("alignment column 4000000")')
    assert b.evaluate('Number(document.querySelector(".alignment-canvas").dataset.firstColumn)>3999900')
    assert b.evaluate('document.querySelector(".alignment-viewport").scrollWidth<=8000000')
    assert b.evaluate('Number(document.querySelector(".alignment-canvas").dataset.paintedCells)<5000')
    close(b)
    # Pagination exposes rows beyond the first DOM window without rendering all cells.
    rows=ROOT/'many.afa';rows.write_text(''.join(f'>r{i}\nAC-GT\n' for i in range(61)))
    open_file(b,rows);fill(b,'select[aria-label="Alignment supplied alphabet"]','Dna','change');b.click('Import alignment');imported(b,4)
    b.evaluate('document.querySelector(".alignment-rows").open=true')
    assert b.evaluate('document.querySelectorAll(".alignment-row").length')==50
    b.click('Next alignment rows');b.wait('document.querySelectorAll(".alignment-row").length===11')
    assert not b.evaluate('document.querySelector(".alignment-panel").textContent.includes("on : click")')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        b.wait('document.querySelector(".alignment-viewport").clientWidth<390')
        assert b.evaluate('document.documentElement.scrollWidth<=390')
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    print('PASS: explicit-alphabet alignment import/column and gap maps, source/export losses, transactional row extraction, native persistence/reopen, stale import rejection, CLUSTAL evidence, 4M-column virtualized scrolling and row pages/mobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))

if __name__=='__main__':main()
