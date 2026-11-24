#!/usr/bin/env python3
"""Actual F25 WASM geometry, selected-reference search and independent donors.
Synthetic sequences and independent literal/IUPAC oracles; no scientific executable.
"""
import base64
import io
import json
import time
import zipfile
from pathlib import Path
from browser_qa import QA_ORIGIN, Browser, mobile_enabled
from methods_browser_qa import record,fill,check,active
from workspace_browser_qa import wait_file,manifest,slots

DIRECTORY=Path('/tmp/collagene2-crispr-qa')
ORIGIN=QA_ORIGIN
MASK={'A':1,'C':2,'G':4,'T':8,'R':5,'Y':10,'S':6,'W':9,'K':12,'M':3,'B':14,'D':13,'H':11,'V':7,'N':15}
COMP=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')


def reverse(s):return s.translate(COMP)[::-1]


def panel(b):b.click('CRISPR');b.wait('!!document.querySelector(".crispr-panel")')


def close(b):
    b.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()')
    b.wait('!document.querySelector("[role=dialog]")')


def paste(b,text):
    b.click('Paste sequence');fill(b,'Sequence input',text);b.click('Open sequence');b.wait('!document.querySelector("[role=dialog]")')


def open_file(b,path,selector='input[type=file]'):
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId']
    node=b.call('DOM.querySelector',{'nodeId':root,'selector':selector})['nodeId'];assert node
    b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})


def archive(b):
    path=DIRECTORY/'collagene-workspace.cg2workspace';path.unlink(missing_ok=True)
    b.wait('document.body?.innerText.includes("Saved in this browser")');b.click('Workspace');b.click('Save workspace file')
    value=wait_file(path);close(b);return value


def evidence(b,label,filename):
    path=DIRECTORY/filename;path.unlink(missing_ok=True);b.click(label);return json.loads(wait_file(path))


def design(b,name,preset='spcas9'):
    active(b,name);panel(b);fill(b,'CRISPR nuclease',preset,'change')
    fill(b,'CRISPR minimum GC','0');fill(b,'CRISPR maximum GC','100');fill(b,'CRISPR homopolymer limit','64');fill(b,'CRISPR candidate limit','4096');fill(b,'CRISPR guide budget','20000000')
    b.click('Design CRISPR guides');b.wait('!!document.querySelector(".crispr-guides")',timeout=60)
    return evidence(b,'Download guide report','crispr-guides.json')


def choose(b,r,spacer,reverse_value=None,starts=None):
    gs=[g for g in r['report']['candidates']if g['guide_sequence']==spacer and g['passes_filters'] and(reverse_value is None or g['reverse']==reverse_value)and(starts is None or g['guide_positions'][0]==starts)]
    assert gs,r
    g=gs[0]
    b.evaluate('(()=>{const e=Array.from(document.querySelectorAll(".crispr-guide-entry")).find(e=>e.dataset.guideId==='+json.dumps(g['id'])+');if(!e)throw Error("Expected guide missing from page");e.click()})()')
    b.wait('document.querySelector(".crispr-guide-detail")?.textContent.includes('+json.dumps(spacer)+')')
    return g


def guide_oracle(r,sequence):
    p=r['options']['profile'];n=len(sequence)
    for g in r['candidates']:
        guide=''.join((sequence[i].translate(COMP)if g['reverse']else sequence[i])for i in g['guide_positions'])
        pam=''.join((sequence[i].translate(COMP)if g['reverse']else sequence[i])for i in g['pam_positions'])
        assert guide==g['guide_sequence'] and pam==g['pam_sequence']
        assert len(g['guide_positions'])==p['guide_length'] and len(g['pam_positions'])==len(p['pam'])
        assert all(MASK[a]&MASK[z]for a,z in zip(pam,p['pam']))
        assert g['concrete']==all(a in 'ACGT'for a in guide)
        assert g['pam_certain']==all(MASK[a]&~MASK[z]==0 for a,z in zip(pam,p['pam']))
        assert all(0<=i<n for i in g['guide_positions']+g['pam_positions'])
        if g['concrete']:assert g['gc_bp']==sum(a in 'GC'for a in guide)*10000//len(guide)


def brute_hits(sequence,circular,profile,query,k):
    n=len(sequence);length=profile['guide_length'];p=len(profile['pam']);window=length+p
    if window>n:return set()
    found=set()
    for rev in [False,True]:
        oriented=reverse(sequence)if rev else sequence
        text=oriented+oriented if circular else oriented
        for start in range(n if circular else n-window+1):
            go=start if profile['pam_side']=='ThreePrime'else start+p
            po=start+length if profile['pam_side']=='ThreePrime'else start
            guide=text[go:go+length];pam=text[po:po+p]
            if not all(MASK[a]&MASK[z]for a,z in zip(pam,profile['pam'])):continue
            minimum=sum(not(MASK[a]&MASK[z])for a,z in zip(guide,query))
            maximum=sum(MASK[a]!=MASK[z]for a,z in zip(guide,query))
            if minimum>k:continue
            indices=lambda at,count:tuple((n-1-(at+i)%n)if rev else(at+i)%n for i in range(count))
            certain=all(MASK[a]&~MASK[z]==0 for a,z in zip(pam,profile['pam']))
            found.add((indices(go,length),indices(po,p),rev,minimum,maximum,certain,maximum<=k and certain))
    return found


def search_refs(b,names):
    b.click('Clear CRISPR references');b.evaluate('document.querySelector(".crispr-reference-picker").open=true')
    for name in names:
        b.evaluate('(()=>{const row=Array.from(document.querySelectorAll(".crispr-reference-choice")).find(e=>e.querySelector("strong").textContent==='+json.dumps(name)+');if(!row)throw Error("Reference missing");row.querySelector("input").click()})()')
    b.wait('document.querySelectorAll(".crispr-reference-label").length==='+str(len(names)))
    labels=b.evaluate('Array.from(document.querySelectorAll(".crispr-reference-label input[aria-label^=\\"CRISPR build\\"]")).map(e=>e.getAttribute("aria-label"))')
    for i,label in enumerate(labels):fill(b,label,'synthetic-selected-record-'+str(i))
    # Updating a label must preserve the focused input rather than recreate it.
    b.evaluate('document.querySelector(".crispr-reference-label input").focus()')
    b.call('Input.insertText',{'text':'-typed'})
    assert b.evaluate('document.activeElement?.getAttribute("aria-label")?.startsWith("CRISPR build")')
    check(b,'CRISPR mark intended locus',True)
    fill(b,'CRISPR maximum mismatches','0');fill(b,'CRISPR search budget','100000000');fill(b,'CRISPR hit limit','4096');fill(b,'CRISPR index limit','1000000')
    b.click('Search CRISPR references');b.wait('!!document.querySelector(".crispr-search-result")',timeout=60)
    return evidence(b,'Download off-target report','crispr-off-targets.json')


def donor_options(b,target,start,deleted,replacement,left,right,distance=100,reversed_donor=False):
    b.evaluate('(()=>{const e=document.querySelector("[aria-label=\\"CRISPR independent target\\"]"),o=Array.from(e.options).find(o=>o.textContent.startsWith('+json.dumps(target+' · ')+'));if(!o)throw Error("Target missing");e.value=o.value;e.dispatchEvent(new Event("change",{bubbles:true}))})()')
    for label,value in [('CRISPR edit start',start+1),('CRISPR delete length',deleted),('CRISPR replacement',replacement),('CRISPR left arm',left),('CRISPR right arm',right),('CRISPR cut distance',distance)]:fill(b,label,str(value))
    check(b,'CRISPR reverse donor',reversed_donor)
    b.click('Design and check donor');b.wait('!!document.querySelector(".crispr-donor-result")',timeout=60)
    return evidence(b,'Download donor review JSON','crispr-donor-review.json')


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for p in DIRECTORY.iterdir():
        if p.is_file():p.unlink()
    b=Browser();b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':ORIGIN});b.wait('document.body?.innerText.includes("Saved in this browser")')
    spacer='ACGTTAGCTAGTCAGTACGA';cas12spacer='ACGTCAGTACGATCGTACGTAGC'
    original='AACCAA'+spacer+'TGG'+'AATTCCAATTCCAATT';edit=23;deleted=2;replacement='TTA';edited=original[:edit]+replacement+original[edit+deleted:]
    circular=original[14:]+original[:14];circular_target='GGA'+circular[2:-2]
    extra=spacer+'TGG'+'AA'+spacer+'TGG';ambiguous=spacer[:5]+'N'+spacer[6:]+'TGG';ambig_pam=spacer+'TGN';many=(spacer+'TGGAAA')*80
    rows=[('cas_source',original,False),('cas_target',edited,False),('cas_wrong',original,False),('extra_ref',extra,False),('ambiguous_ref',ambiguous,False),('ambiguous_pam',ambig_pam,False),('cas12_source','AACCAA'+'TTTA'+cas12spacer+'ACAACACACA',False),('reverse_source',reverse(original),False),('circle_source',circular,True),('circle_target',circular_target,True),('many_reference',many,False),('custom_source','AACCAA'+spacer+'GG'+'AATTAA',False),('poly_source','A'*20+'TGG',False),('cas_source',original,False)]
    sequences={name:(seq,circ)for name,seq,circ in rows}
    paste(b,''.join(record(*row)for row in rows))
    r=design(b,'cas_source');guide_oracle(r['report'],original);g=choose(b,r,spacer,False,6)
    assert g['pam_sequence']=='TGG' and g['nominal_top_cut']==23
    refs=search_refs(b,['cas_source','extra_ref','ambiguous_ref','ambiguous_pam'])
    sr=refs['report'];assert sr['complete'] and sr['exhaustive_for_selected_records'];assert len(refs['reference_events'])==4
    for i,reference in enumerate(sr['references']):
        seq,circ=sequences[reference['name']]
        actual={(tuple(h['guide_positions']),tuple(h['pam_positions']),h['reverse'],h['min_mismatches'],h['max_mismatches'],h['pam_certain'],h['definite'])for h in sr['hits']if h['reference_index']==i}
        assert actual==brute_hits(seq,circ,sr['options']['profile'],spacer,0),(reference,actual)
    assert sum(h['is_intended']for h in sr['hits'])==1
    assert sum(h['min_mismatches']==h['max_mismatches']==0 and not h['is_intended']for h in sr['hits'])>=2
    assert any(h['min_mismatches']==0 and h['max_mismatches']==1 and not h['definite']for h in sr['hits'])
    assert any(not h['pam_certain']for h in sr['hits'])
    b.evaluate('document.querySelector(".crispr-hit-entry").click()');b.wait('!!document.querySelector(".crispr-hit-detail")')
    b.click('Open reference guide position');b.wait('!document.querySelector(".crispr-panel")');b.wait('document.querySelector(".statusbar")?.textContent.includes("Selection 7–7")')
    # Both nominal profiles, reverse coordinates and original-circle wrap mappings.
    cr=design(b,'cas12_source','cas12a');guide_oracle(cr['report'],sequences['cas12_source'][0]);cg=choose(b,cr,cas12spacer,False,10);assert cg['pam_sequence']=='TTTA'
    rr=design(b,'reverse_source');guide_oracle(rr['report'],reverse(original));rg=choose(b,rr,spacer,True);assert rg['guide_positions'][0]>rg['guide_positions'][-1]
    cc=design(b,'circle_source');guide_oracle(cc['report'],circular);choose(b,cc,spacer,False);assert any(a>b for a,b in zip(next(g for g in cc['report']['candidates']if g['guide_sequence']==spacer and not g['reverse'])['guide_positions'],next(g for g in cc['report']['candidates']if g['guide_sequence']==spacer and not g['reverse'])['guide_positions'][1:]))
    cir=donor_options(b,'circle_target',len(circular)-2,4,'GGA',5,5,20)
    expected_circle=circular[-7:-2]+'GGA'+circular[2:7]
    assert cir['independent_check']['consistent'] and cir['report']['donor']['sequence']['bases']==expected_circle
    assert cir['report']['edited_reference']['sequence']['bases']==circular_target
    # Explicit custom geometry remains declared, without borrowing preset authority.
    cu=design(b,'custom_source','custom')
    # Custom default PAM NGG does not match the desired GG-only locus. Change it explicitly.
    fill(b,'CRISPR custom PAM','GG');b.click('Design CRISPR guides');b.wait('!!document.querySelector(".crispr-guides")')
    cu=evidence(b,'Download guide report','crispr-guides.json');guide_oracle(cu['report'],sequences['custom_source'][0]);choose(b,cu,spacer,False,6);assert cu['report']['options']['profile']['custom']
    poly=design(b,'poly_source');fill(b,'CRISPR homopolymer limit','3');b.click('Design CRISPR guides');b.wait('!!document.querySelector(".crispr-guides")')
    pr=evidence(b,'Download guide report','crispr-guides.json');assert any(not g['passes_filters']for g in pr['report']['candidates'])
    check(b,'CRISPR show rejected guides');assert b.evaluate('document.querySelectorAll(".crispr-guide-entry").length')>0
    # Candidate/hit limits never establish absence; retained viewers remain bounded.
    mr=design(b,'many_reference');assert len(mr['report']['candidates'])>25
    assert b.evaluate('document.querySelectorAll(".crispr-guide-entry").length')==25
    b.click('Next guides');assert b.evaluate('document.querySelectorAll(".crispr-guide-entry").length')<=25
    fill(b,'CRISPR candidate limit','1');b.click('Design CRISPR guides');b.wait('!!document.querySelector(".crispr-guides")');limited=evidence(b,'Download guide report','crispr-guides.json');assert not limited['report']['complete'] and limited['report']['limits']
    r=design(b,'cas_source');choose(b,r,spacer,False,6)
    hits=search_refs(b,['cas_source','many_reference']);assert hits['report']['complete'] and len(hits['report']['hits'])>=81
    assert b.evaluate('document.querySelectorAll(".crispr-hit-entry").length')==25
    b.click('Next off-target hits');assert b.evaluate('document.querySelectorAll(".crispr-hit-entry").length')==25
    fill(b,'CRISPR hit limit','1');b.click('Search CRISPR references');b.wait('!!document.querySelector(".crispr-search-result")');lim=evidence(b,'Download off-target report','crispr-off-targets.json');assert not lim['report']['complete'] and lim['report']['limits']
    fill(b,'CRISPR hit limit','4096');fill(b,'CRISPR search budget','1');b.click('Search CRISPR references');b.wait('!!document.querySelector(".crispr-search-result")');lim=evidence(b,'Download off-target report','crispr-off-targets.json');assert not lim['report']['complete']
    # A separately supplied wrong target blocks release; literal donor sequence
    # comes from source arms/replacement, never from that target.
    wrong=donor_options(b,'cas_wrong',edit,deleted,replacement,6,6,10);assert not wrong['independent_check']['consistent'] and not wrong['independent_check']['target_matches']
    assert b.evaluate('Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Open checked donor template").disabled')
    good=donor_options(b,'cas_target',edit,deleted,replacement,6,6,10);expected=original[edit-6:edit]+replacement+original[edit+deleted:edit+deleted+6]
    assert good['independent_check']['consistent'] and good['report']['donor']['sequence']['bases']==expected
    assert good['report']['edited_reference']['sequence']['bases']==edited
    rev=donor_options(b,'cas_target',edit,deleted,replacement,6,6,10,True);assert rev['report']['donor']['sequence']['bases']==reverse(expected) and rev['independent_check']['consistent']
    good=donor_options(b,'cas_target',edit,deleted,replacement,6,6,10)
    b.click('Export checked donor FASTA');fasta=wait_file(DIRECTORY/'crispr-donor.fasta').decode();assert ''.join(fasta.splitlines()[1:])==expected
    b.click('Export checked donor GenBank');gb=wait_file(DIRECTORY/'crispr-donor.gb').decode();assert 'ORIGIN' in gb and 'FEATURES' in gb
    before_data=archive(b);before=manifest(before_data);before_slots=slots(b)
    # Atomic second-put failure cannot install half a donor or mutate source heads.
    b.click('Open checked donor template');b.wait('!!document.querySelector(".crispr-open-panel")',timeout=60)
    b.evaluate("window.qaPut=IDBObjectStore.prototype.put;IDBObjectStore.prototype.put=function(v,k){if(k==='current')throw new DOMException('Synthetic donor transaction failure','QuotaExceededError');return window.qaPut.call(this,v,k)};true")
    b.click('Recheck and open CRISPR product');b.wait('document.querySelector(".crispr-open-panel")?.textContent.includes("Donor product was not opened")',timeout=60)
    b.evaluate('IDBObjectStore.prototype.put=window.qaPut;delete window.qaPut;true');assert slots(b)==before_slots
    b.click('Recheck and open CRISPR product');b.wait('!document.querySelector(".crispr-open-panel")',timeout=60)
    saved_data=archive(b);saved=manifest(saved_data);assert len(saved['documents'])==len(before['documents'])+1
    product=next(d for d in saved['documents']if d['id']==saved['active_document']);event=next(e for e in saved['history']if e['id']==product['head_event'])
    assert event['parents']==[good['source']['event_id']]
    assert good['independent_target']['event_id'] not in event['parents']
    for old in before['documents']:assert next(d for d in saved['documents']if d['id']==old['id'])==old
    assert manifest(base64.b64decode(slots(b)['previous']))['documents']==before['documents']
    assert 'unverified on opening' in b.evaluate('document.querySelector(".crispr-historical").textContent')
    b.click('Download retained CRISPR provenance');retained=json.loads(wait_file(DIRECTORY/'crispr-retained-provenance.json'));assert retained['source']==good['source'] and retained['independent_target']==good['independent_target']
    portable=DIRECTORY/'crispr-native.cg2workspace';portable.write_bytes(saved_data)
    b.click('Workspace');open_file(b,portable,'[role=dialog] input[type=file][accept*=cg2workspace]');b.wait('!document.querySelector("[role=dialog]")',timeout=60)
    b.call('Page.reload');b.wait('!!document.querySelector(".crispr-historical")')
    assert b.evaluate('document.querySelector(".crispr-donor-result")===null')
    # Every checked edited reference also has only source ancestry; target remains
    # a noncausal requirement and current verification is rebuilt from scratch.
    r=design(b,'cas_source');choose(b,r,spacer,False,6);good=donor_options(b,'cas_target',edit,deleted,replacement,6,6,10)
    b.click('Open checked edited reference');b.wait('!!document.querySelector(".crispr-open-panel")',timeout=60);b.click('Recheck and open CRISPR product');b.wait('!document.querySelector(".crispr-open-panel")',timeout=60)
    ed=manifest(archive(b));product=next(d for d in ed['documents']if d['id']==ed['active_document']);assert next(e for e in ed['history']if e['id']==product['head_event'])['parents']==[good['source']['event_id']]
    # Synchronous options, equal-content document identity, and explicit cancel
    # guards terminate pending workers before replies can become current.
    r=design(b,'cas_source');choose(b,r,spacer,False,6)
    js='''(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Design CRISPR guides").click();const e=document.querySelector('[aria-label="CRISPR minimum GC"]');e.value="10";e.dispatchEvent(new Event("input",{bubbles:true}));})()'''
    b.evaluate(js);time.sleep(.3);assert b.evaluate('document.querySelector(".crispr-guides")===null')
    r=design(b,'cas_source')
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Design CRISPR guides").click();queueMicrotask(()=>Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Cancel guide design").click())})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".crispr-guides")===null')
    r=design(b,'cas_source')
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Design CRISPR guides").click();Array.from(document.querySelectorAll(".document-entry")).filter(e=>e.querySelector("strong").textContent==="cas_source").at(-1).click()})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".crispr-guides")===null')
    r=design(b,'cas_source');choose(b,r,spacer,False,6);donor_options(b,'cas_target',edit,deleted,replacement,6,6,10)
    b.click('Expand CRISPR');b.wait('document.querySelector(".crispr-panel").classList.contains("crispr-expanded")')
    b.evaluate('document.querySelector(".crispr-donor-result").scrollIntoView({block:"center"})');b.screenshot('/tmp/collagene2-crispr-desktop.png')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        b.evaluate('document.querySelector(".crispr-donor-result").scrollIntoView({block:"center"})');assert b.evaluate('document.documentElement.scrollWidth<=innerWidth');b.screenshot('/tmp/collagene2-crispr-mobile.png')
    assert not b.exceptions and not b.console_errors,(b.exceptions,b.console_errors)
    print('PASS: SpCas9/Cas12a/custom geometry, independent circular/reverse/IUPAC coordinates and exact selected-record Hamming oracle, intended/additional matches, limits/paging/cancel, literal wrong-target guarded donors, source-only parents, atomic/native provenance and mobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__=='__main__':main()
