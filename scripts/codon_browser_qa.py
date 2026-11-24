#!/usr/bin/env python3
"""Actual WASM F26 design/check workflow against independent literal peptides.
Synthetic sequences only; run with isolated Chrome :9224 and preview :8080.
"""
import base64
import json
import time
from pathlib import Path

from browser_qa import Browser, mobile_enabled
from methods_browser_qa import fill, record, check
from workspace_browser_qa import slots, wait_file, manifest

DIRECTORY = Path('/tmp/collagene2-codon-qa')
ORIGIN = 'http://127.0.0.1:8080'


def activate(b, name, occurrence=0):
    b.evaluate('Array.from(document.querySelectorAll(".document-entry")).filter(e=>e.querySelector("strong")?.textContent==='+json.dumps(name)+')['+str(occurrence)+'].click()')
    b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name))


def protein(b, name, sequence):
    b.click('Paste sequence'); fill(b,'Sequence input',f'>{name} [molecule=protein]\n{sequence}')
    b.click('Open sequence')
    title=name+' [molecule=protein]'
    b.wait('document.querySelector("h1")?.textContent==='+json.dumps(title))
    b.click('Codon design'); b.wait('!!document.querySelector(".codon-design-panel")')
    return title


def design(b, status='CandidateFound'):
    b.click('Design coding DNA')
    b.wait('document.querySelector(".codon-design-result")?.dataset.status==='+json.dumps(status),timeout=60)
    if status=='CandidateFound':
        b.wait('document.querySelector(".codon-check-result")?.dataset.satisfies==="true"',timeout=60)


def report(b):
    path=DIRECTORY/'codon-design-report.json';path.unlink(missing_ok=True)
    b.click('Download codon report');return json.loads(wait_file(path))


def archive(b):
    path=DIRECTORY/'collagene-workspace.cg2workspace';path.unlink(missing_ok=True)
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.click('Workspace');b.wait('!!document.querySelector("[role=dialog]")')
    b.click('Save workspace file');data=wait_file(path)
    b.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()')
    b.wait('!document.querySelector("[role=dialog]")');return data


def check_supplied(b, dna, passing):
    fill(b,'Codon supplied DNA',dna)
    b.click('Check coding DNA')
    b.wait('!!document.querySelector(".codon-check-result")',timeout=60)
    assert b.evaluate('document.querySelector(".codon-check-result").dataset.satisfies')==str(passing).lower()
    assert b.evaluate('Array.from(document.querySelectorAll("button")).some(e=>e.textContent==="Open designed DNA")')==passing


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for path in DIRECTORY.iterdir():
        if path.is_file():path.unlink()
    b=Browser()
    b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':ORIGIN})
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    title=protein(b,'literal_codon_source','MKW*')
    fill(b,'Codon weights','AAG 100')
    design(b)
    r=report(b);d=r['design'];c=r['independent_check']
    assert d['candidate']['sequence']['bases']=='ATGAAGTGGTAA',d
    assert d['target_protein']==c['target_protein']==c['translated_protein']=='MKW*'
    assert c['score']==103 and d['search_complete'] and d['optimality_proven']
    assert b.evaluate('document.querySelector("h1").textContent')==title,'Design silently replaced its source'
    assert b.evaluate('document.querySelector(".codon-check-result").dataset.complete')=='true'
    # A forbidden non-palindrome also excludes its reverse complement AAG.
    fill(b,'Codon forbidden motifs','CTT')
    assert b.evaluate('document.querySelector(".codon-design-result")===null')
    design(b);r=report(b);assert r['design']['candidate']['sequence']['bases']=='ATGAAATGGTAA'
    # The unique initial ATG makes this junction constraint impossible.
    fill(b,'Codon forbidden motifs','CATG');fill(b,'Codon left context','C')
    design(b,'Unsatisfiable');r=report(b)
    assert r['design']['search_complete'] and r['design']['candidate'] is None
    assert not b.evaluate('Array.from(document.querySelectorAll("button")).some(e=>e.textContent==="Open designed DNA")')
    fill(b,'Codon forbidden motifs','');fill(b,'Codon left context','')
    fill(b,'Codon minimum GC','100');design(b,'Unsatisfiable')
    fill(b,'Codon minimum GC','0');fill(b,'Codon weights','AAG 0');fill(b,'Codon maximum homopolymer','2')
    design(b,'Unsatisfiable')
    fill(b,'Codon weights','AAG 100');fill(b,'Codon maximum homopolymer','12')
    check(b,'Codon window GC enabled',True);fill(b,'Codon GC window','1');fill(b,'Codon window maximum GC','0')
    design(b,'Unsatisfiable')
    check(b,'Codon window GC enabled',False);fill(b,'Codon work limit','1')
    design(b,'SearchIncomplete');r=report(b)
    assert not r['design']['search_complete'] and r['design']['candidate'] is None
    fill(b,'Codon work limit','2000000')
    # The checker is available independently of the optimizer and rejects a
    # one-codon amino-acid change even though the DNA is syntactically valid.
    check_supplied(b,'ATGAAGTGGTAA',True)
    check_supplied(b,'ATGGAGTGGTAA',False)
    r=report(b);assert r['independent_check']['translated_protein']=='MEW*'
    fill(b,'Codon supplied DNA','')
    # Concrete DNA sources preserve internal-codon meanings; code11 GTG is V
    # here rather than silently asserting bacterial translation initiation.
    b.click('Paste sequence');fill(b,'Sequence input',record('ordinary_valine','GTG'))
    b.click('Open sequence');b.wait('document.querySelector("h1")?.textContent==="ordinary_valine"')
    fill(b,'Codon genetic code','11','change');fill(b,'Codon weights','GTA 100')
    design(b);r=report(b)
    assert r['design']['target_protein']=='V' and r['design']['candidate']['sequence']['bases']=='GTA'
    # Mitochondrial ATA and TGA are independently specified M and W.
    mito=protein(b,'mitochondrial_codon_source','MW')
    fill(b,'Codon genetic code','2','change');fill(b,'Codon weights','ATA 100\nTGA 100')
    design(b);r=report(b)
    assert r['design']['candidate']['sequence']['bases']=='ATATGA' and r['independent_check']['translated_protein']=='MW'
    fill(b,'Codon genetic code','1','change');fill(b,'Codon weights','')
    long_title=protein(b,'bounded_leucine_source','LLLLLLLLLL')
    fill(b,'Codon beam width','1');design(b)
    r=report(b);assert r['design']['pruned_states']>0 and not r['design']['optimality_proven'] and not r['design']['search_complete']
    # Immediately changing constraints invalidates an asynchronous request.
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Design coding DNA").click();const e=document.querySelector("[aria-label=\\"Codon maximum homopolymer\\"]");e.value="10";e.dispatchEvent(new Event("input",{bubbles:true}));})()')
    time.sleep(.3)
    assert b.evaluate('document.querySelector(".codon-check-result")===null')
    assert not b.evaluate('Array.from(document.querySelectorAll("button")).some(e=>e.textContent==="Open designed DNA")')
    # Duplicate identical-content documents cannot inherit a pending result.
    duplicate=protein(b,'literal_codon_source','MKW*')
    fill(b,'Codon beam width','128');fill(b,'Codon maximum homopolymer','12');fill(b,'Codon weights','AAG 100')
    b.evaluate('Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Design coding DNA").click()')
    activate(b,title,0);time.sleep(.3)
    assert b.evaluate('document.querySelector(".codon-check-result")===null')
    design(b)
    before=manifest(archive(b));source_doc=next(d for d in before['documents'] if d['id']==before['active_document'])
    b.click('Download designed GenBank')
    genbank=wait_file(DIRECTORY/'codon-designed-dna.gb').decode()
    assert 'ATGAAGTGGTAA'.lower() in ''.join(genbank.split()).lower()
    assert 'codon' in genbank.lower() and 'MKW*' in genbank
    b.screenshot('/tmp/collagene2-codon-desktop.png')
    b.click('Open designed DNA')
    b.wait('document.querySelector("h1")?.textContent!=='+json.dumps(title))
    after=manifest(archive(b))
    assert next(d for d in after['documents'] if d['id']==source_doc['id'])==source_doc
    product=next(d for d in after['documents'] if d['id']==after['active_document'])
    event=next(e for e in after['history'] if e['id']==product['head_event'])
    assert event['parents']==[source_doc['head_event']],event
    # Native recovery retains provenance, but never restores a current report.
    b.call('Page.reload');b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.click('Codon design');assert b.evaluate('document.querySelector(".codon-check-result")===null')
    slots_after=slots(b);assert slots_after['current']
    saved=manifest(base64.b64decode(slots_after['current']))
    assert any(e['id']==event['id'] and e['parents']==event['parents'] for e in saved['history'])
    # Reopening the exported GenBank preserves declared coding evidence.
    path=DIRECTORY/'codon-designed-dna.gb'
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId']
    node=b.call('DOM.querySelector',{'nodeId':root,'selector':'input[type=file][multiple]'})['nodeId']
    b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.wait('document.querySelectorAll(".document-entry").length==='+str(len(after['documents'])+1))
    activate(b,title,0);b.click('Codon design');design(b)
    full_gene=protein(b,'full_gene_codon_source','L'*1000)
    fill(b,'Codon beam width','128');fill(b,'Codon work limit','20000000');fill(b,'Codon weights','')
    design(b);long_report=report(b)
    assert len(long_report['design']['candidate']['sequence']['bases'])==3000
    assert long_report['independent_check']['translated_protein']=='L'*1000
    assert long_report['independent_check']['satisfies'] is True and long_report['independent_check']['complete']
    assert b.evaluate('document.querySelector("h1").textContent')==full_gene,'Full-gene design replaced source'
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        assert b.evaluate('document.documentElement.scrollWidth<=innerWidth'),'Codon panel overflows mobile viewport'
        b.screenshot('/tmp/collagene2-codon-mobile.png')
    b.click('Design coding DNA');b.click('Features');time.sleep(.3)
    assert not b.exceptions and not b.console_errors,(b.exceptions,b.console_errors)
    print('PASS: literal codon preferences/translation, both-strand/context constraints, exhaustive vs incomplete search, independent DNA checks, alternate code/ordinary initiation, pruning/stale-source cancellation, exact derived parent, native/GenBank provenance and mobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__=='__main__':main()
