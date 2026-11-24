#!/usr/bin/env python3
"""Actual checked positional restriction/repair/ligation, retention and replay QA.
Only literal synthetic substrates/targets. External pdftotext is a dev export oracle.
"""
import io
import json
import subprocess
import time
import zipfile
from pathlib import Path
from browser_qa import Browser, mobile_enabled
from methods_browser_qa import record,fill,check,source,active
from workspace_browser_qa import wait_file,manifest

DIRECTORY=Path('/tmp/collagene2-checked-restriction-qa')
ORIGIN='http://127.0.0.1:8080'


def close(b):
    b.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()')
    b.wait('!document.querySelector("[role=dialog]")')


def archive(b):
    path=DIRECTORY/'collagene-workspace.cg2workspace';path.unlink(missing_ok=True)
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.click('Workspace');b.click('Save workspace file');data=wait_file(path);close(b);return data


def open_file(b,path,selector='input[type=file]'):
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId']
    node=b.call('DOM.querySelector',{'nodeId':root,'selector':selector})['nodeId'];assert node
    b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})


def run_records(data):
    m=manifest(data)
    with zipfile.ZipFile(io.BytesIO(data))as z:
        return [json.loads(z.read('blobs/'+a['blob']))for a in m['attachments'].values()if a['kind']=='cloning-run-v1']


def report(b):
    path=DIRECTORY/'cloning-reconstruction.json';path.unlink(missing_ok=True)
    b.click('Download reconstruction report');return json.loads(wait_file(path))


def states(b,dam='unmodified',dcm='unmodified',cpg='unmodified',side='First'):
    for category,value in [('Dam',dam),('Dcm',dcm),('CpG',cpg)]:fill(b,f'{side} checked {category} assumption',value,'change')


def configure(b,target,substrate,enzyme,dam='unmodified',dcm='unmodified',cpg='unmodified',repair='none',policy='checked',final='ligate'):
    active(b,target);b.click('Plan')
    if b.evaluate('!!document.querySelector(".candidate-review")'):b.click('Discard candidate')
    source(b,'First source document',substrate);fill(b,'Cloning method','restriction','change')
    fill(b,'Restriction model',policy,'change');fill(b,'Restriction final product',final,'change')
    if final=='ligate':check(b,'Use second restriction substrate',False)
    fill(b,'First enzymes',enzyme);fill(b,'First fragment number','1');fill(b,'First end repair',repair,'change');fill(b,'First end treatment','none','change');check(b,'Reverse First fragment',False)
    if final=='ligate':check(b,'Allow restriction ligation nicks',False)
    if policy=='checked':states(b,dam,dcm,cpg)


def reconstruct(b):
    b.click('Reconstruct and verify');b.wait('!!document.querySelector(".verification-result")',timeout=60)
    return report(b)


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for p in DIRECTORY.iterdir():
        if p.is_file():p.unlink()
    b=Browser();b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':ORIGIN});b.wait('document.body?.innerText.includes("Saved in this browser")')
    rows=[('dam_source','TCTAGATCAAA',True),('dam_target','TCTAGATCAAA',True),('wrong_target','CCCGGG',True),('fill_source','GAATTCAAAA',True),('fill_target','AATTCAAAAGAATT',True),('cpg_source','CCCGGGAAAA',True),('cpg_target','CCCGGGAAAA',True),('tracked_fill_source','GAATTCAAAAGATC',True),('tracked_fill_target','AATTCAAAAGATCGAATT',True),('dam_source','TCTAGATCAAA',True),('dam_target','TCTAGATCAAA',True)]
    b.click('Paste sequence');fill(b,'Sequence input',''.join(record(*r)for r in rows));b.click('Open sequence')
    b.wait('!document.querySelector("[role=dialog]")')
    # One checked, fully Dam-modified XbaI circle stays uncut and cannot become
    # a circular ligation substrate with invented free ends.
    configure(b,'dam_target','dam_source','XbaI',dam='fully',final='fragment')
    unchanged=reconstruct(b);assert unchanged['status']=='ConsistentUnderAssumptions' and len(unchanged['products'])==2
    assert unchanged['products'][-1]['topology']=='Circular' and unchanged['products'][-1]['sequence']['bases']=='TCTAGATCAAA'
    assert 'blocked' in ' '.join(unchanged['findings']).lower()
    fill(b,'Restriction final product','ligate','change');check(b,'Use second restriction substrate',False)
    blocked=reconstruct(b);assert blocked['status']!='ConsistentUnderAssumptions'
    text=' '.join(blocked['findings']).lower();assert 'blocked'in text
    assert len(blocked['products'])<3,'Blocked circle must not gain a final ligated product'
    # Unknown eligibility withholds the digest rather than selecting an arbitrary
    # subset of possible products. Incomplete-run acknowledgement cannot invent it.
    states(b,dam='unknown');uncertain=reconstruct(b);assert uncertain['status']=='Incomplete'
    assert len(uncertain['products'])==1
    b.click('Keep reconstructed products');b.wait('!!document.querySelector(".run-capture-panel")')
    check(b,'Allow incomplete run recording',True);b.click('Record cloning run')
    b.wait('document.querySelector(".run-capture-panel")?.textContent.includes("Run was not recorded")',timeout=60)
    assert 'complete reconstruction' in b.evaluate('document.querySelector(".run-capture-panel").textContent')
    close(b)
    states(b,dam='mixed');mixed=reconstruct(b);assert mixed['status']=='Incomplete'
    # A literal unmodified circle can be cut and self-ligated; comparison remains
    # against the separately selected independent target document.
    states(b);permitted=reconstruct(b);assert permitted['status']=='ConsistentUnderAssumptions'
    assert len(permitted['products'])==3
    sequence=permitted['products'][-1]['sequence']['bases'];assert sequence in ('TCTAGATCAAA'*2),'Circular gauge must be an observed rotation'
    assert len(sequence)==11
    for product in permitted['products']:
        summary=json.loads(product['metadata']['tracked_modifications'])
        assert summary['known_modified']==0 and summary['unknown_eligible']==0
    # Equal-content source and target documents retain distinct event identities.
    # Their pending results cannot cross a source choice or active-target switch.
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Reconstruct and verify").click();const e=document.querySelector(\'[aria-label="First source document"]\');e.value=Array.from(e.options).filter(o=>o.textContent==="dam_source").at(-1).value;e.dispatchEvent(new Event("change",{bubbles:true}));})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".verification-result")===null')
    configure(b,'dam_target','dam_source','XbaI');reconstruct(b)
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Reconstruct and verify").click();Array.from(document.querySelectorAll(".document-entry")).filter(e=>e.querySelector("strong").textContent==="dam_target").at(-1).click();})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".verification-result")===null')
    configure(b,'dam_target','dam_source','XbaI');permitted=reconstruct(b)
    before_data=archive(b);before=manifest(before_data)
    b.click('Keep reconstructed products');b.wait('!!document.querySelector(".run-capture-panel")')
    selected=b.evaluate('document.querySelector("[aria-label=\\"Source event for first\\"]").value');assert selected
    source_doc=next(d for d in before['documents']if d['id']==next(e for e in before['history']if e['id']==selected)['document_id'])
    assert source_doc['head_event']==selected
    b.click('Record cloning run');b.wait('!document.querySelector(".run-capture-panel")',timeout=60)
    saved_data=archive(b);saved=manifest(saved_data);records=run_records(saved_data);assert len(records)==1
    run=records[0];assert [next(iter(s))for s in run['plan']['steps']]==['DeclareSourceModifications','CheckedDigest','Ligate']
    assert len(run['inputs'])==1 and run['inputs'][0]['event']['event_id']==selected
    assert saved['active_document']==before['active_document'],'Independent target remains active'
    events={e['id']:e for e in saved['history']}
    assert events[run['products'][0]['event_id']]['parents']==[selected]
    assert events[run['products'][1]['event_id']]['parents']==[run['products'][0]['event_id']]
    assert events[run['products'][2]['event_id']]['parents']==[run['products'][1]['event_id']]
    for d in before['documents']:assert next(n for n in saved['documents']if n['id']==d['id'])==d
    assert saved['plan_schema']==4
    # Native reopening retains historical evidence and the exact operation DAG,
    # then reconstructs positional state on explicit replay.
    portable=DIRECTORY/'checked-native.cg2workspace';portable.write_bytes(saved_data)
    b.click('Workspace');open_file(b,portable,'[role=dialog] input[type=file][accept*=cg2workspace]')
    b.wait('!document.querySelector("[role=dialog]")',timeout=60)
    b.click('History');b.click('View cloning run');b.wait('!!document.querySelector(".run-details-panel")')
    assert 'unverified on opening' in b.evaluate('document.querySelector(".run-authority").textContent')
    b.click('Replay cloning run');b.wait('!!document.querySelector(".run-replay-result")',timeout=60)
    assert 'ConsistentUnderAssumptions' in b.evaluate('document.querySelector(".run-replay-result").textContent')
    b.click('Download run replay JSON');replay=json.loads(wait_file(DIRECTORY/'cloning-run-replay.json'))
    assert replay['products_match'] and replay['status_matches']
    close(b)
    # Checked state continues through literal polymerase fill-in and covalent
    # circular ligation. Newly synthesized calls are unmodified, not methylated
    # merely because an original-generation assumption existed.
    configure(b,'fill_target','fill_source','EcoRI',repair='blunt')
    repaired=reconstruct(b);assert repaired['status']=='ConsistentUnderAssumptions'
    assert len(repaired['products'])==4
    product=repaired['products'][-1];seq=product['sequence']['bases'];assert len(seq)==14 and seq in 'AATTCAAAAGAATT'*2
    mod=json.loads(product['metadata']['tracked_modifications']);assert mod['known_modified']==0 and mod['unknown_eligible']==0
    assert b.evaluate('document.querySelectorAll(".checked-modification-summary").length')==4
    b.click('Export HTML protocol');html=wait_file(DIRECTORY/'cloning-protocol.html').decode()
    assert 'Checked' in html and ('generation'in html.lower()or'modification'in html.lower()) and 'repair'in html.lower()
    b.click('Export PDF protocol');pdf=wait_file(DIRECTORY/'cloning-protocol.pdf');assert pdf.startswith(b'%PDF-')
    pdftext=subprocess.run(['pdftotext',str(DIRECTORY/'cloning-protocol.pdf'),'-'],check=True,capture_output=True,text=True).stdout
    assert 'Checked digest' in pdftext and 'End repair' in pdftext and 'Dam Unmodified' in pdftext
    assert '14 symbols' in pdftext and '0 retained nicks' in pdftext
    assert 'sensitivity profile' in pdftext and 'unknown' in pdftext.lower()
    # Existing Dam calls survive actual-strand fill and ligation; synthesized
    # bases remain unmodified rather than inheriting a source-wide methylase flag.
    configure(b,'tracked_fill_target','tracked_fill_source','EcoRI',dam='fully',repair='blunt')
    tracked=reconstruct(b);assert tracked['status']=='ConsistentUnderAssumptions'
    summaries=[json.loads(p['metadata']['tracked_modifications'])for p in tracked['products']]
    assert all(s['known_modified']==2 and s['unknown_eligible']==0 for s in summaries)
    assert summaries[-1]['top_length']==18 and summaries[-1]['bottom_length']==18
    assert summaries[-1]['known_unmodified']-summaries[0]['known_unmodified']==4
    assert len(tracked['products'][-1]['sequence']['bases'])==18 and tracked['products'][-1]['sequence']['bases'] in 'AATTCAAAAGATCGAATT'*2
    b.click('Save candidate');candidate=wait_file(DIRECTORY/'candidate.cg2plan');assert candidate[4:6]==bytes([6,0])
    # Readable new operation descriptions survive candidate transfer; the target
    # is never silently replaced with the reconstructed output.
    open_file(b,DIRECTORY/'candidate.cg2plan');b.wait('!!document.querySelector(".candidate-review")')
    assert 'Declare original source modifications'in b.evaluate('document.querySelector(".candidate-review").textContent')
    active(b,'wrong_target');b.click('Verify candidate');b.wait('document.querySelector(".notice")?.textContent.includes("active target differs")')
    configure(b,'cpg_target','cpg_source','SmaI',cpg='fully');cp=reconstruct(b);assert cp['status']!='ConsistentUnderAssumptions' and 'blocked' in ' '.join(cp['findings']).lower()
    # An imported producer summary cannot authorize fresh checked state. Fixture
    # is created by root's independent native audit, not by a scientific runtime.
    forged=Path('/tmp/collagene2-checked-forged.cg2');assert forged.exists()
    open_file(b,forged);b.wait('document.querySelector("h1")?.textContent==="forged_checked_source"')
    configure(b,'dam_target','forged_checked_source','XbaI',dam='fully');fr=reconstruct(b)
    assert fr['status']!='ConsistentUnderAssumptions' and 'blocked'in' '.join(fr['findings']).lower()
    # Legacy choice remains explicit and retains its original untracked-cleavage
    # behavior rather than inheriting the checked declaration's state.
    configure(b,'dam_target','dam_source','XbaI',policy='legacy');legacy=reconstruct(b);assert legacy['status']=='ConsistentUnderAssumptions'
    configure(b,'fill_target','fill_source','EcoRI',repair='blunt');reconstruct(b)
    b.evaluate('document.querySelector(".checked-source-assumptions").scrollIntoView({block:"center"})')
    b.screenshot('/tmp/collagene2-checked-restriction-desktop.png')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        b.evaluate('document.querySelector(".checked-source-assumptions").scrollIntoView({block:"center"})')
        assert b.evaluate('document.documentElement.scrollWidth<=innerWidth')
        assert b.evaluate('getComputedStyle(document.querySelector(".checked-restriction-options .checkbox-label")).flexDirection==="row"')
        b.screenshot('/tmp/collagene2-checked-restriction-mobile.png')
    assert not b.exceptions and not b.console_errors,(b.exceptions,b.console_errors)
    print('PASS: literal checked blocked/permitted circles, unknown/unsupported withholding and retention gates, exact source-event parents, positional repair/ligation, native replay, independent target, forged-summary nonauthority, protocols and mobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__=='__main__':main()
