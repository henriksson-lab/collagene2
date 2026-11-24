#!/usr/bin/env python3
"""Actual F21 WASM review, exact binding and bounded site evidence. Synthetic DNA only.
Requires isolated local preview8080 / Chrome9224. No scientific runtime dependency.
"""
import json
import time
from pathlib import Path
from browser_qa import Browser, mobile_enabled
from methods_browser_qa import fill, check, record, active
from workspace_browser_qa import wait_file, manifest

DIRECTORY=Path('/tmp/collagene2-methylation-qa')
ORIGIN='http://127.0.0.1:8080'
MASK={'A':1,'C':2,'G':4,'T':8,'R':5,'Y':10,'S':6,'W':9,'K':12,'M':3,'B':14,'D':13,'H':11,'V':7,'N':15}
COMPLEMENT=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')


def reverse(text):return text.translate(COMPLEMENT)[::-1]


def paste(b,text):
    b.click('Paste sequence');fill(b,'Sequence input',text);b.click('Open sequence')
    b.wait('!document.querySelector("[role=dialog]")')


def panel(b):
    b.click('Methylation');b.wait('!!document.querySelector(".methylation-panel")')


def select(b,names):
    b.click('Clear methylation profiles')
    for name in names:check(b,'Methylation enzyme '+name)


def states(b,dam='unmodified',dcm='unmodified',cpg='unmodified'):
    for label,value in [('Dam',dam),('Dcm',dcm),('CpG',cpg)]:fill(b,'Methylation '+label+' assumption',value,'change')
    fill(b,'Methylation comparison budget','20000000');fill(b,'Methylation site limit','4096')


def assess(b):
    b.click('Assess methylation');b.wait('!!document.querySelector(".methylation-result")',timeout=60)


def report(b):
    path=DIRECTORY/'methylation-assessment.json';path.unlink(missing_ok=True)
    b.click('Download methylation report');value=json.loads(wait_file(path))
    assert value['kind']=='methylation-assessment-review-v1'
    assert value['source_revision']==value['report']['source_revision']
    assert value['source_document_id'] and value['source_event_id']
    return value


def close(b):
    b.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()')
    b.wait('!document.querySelector("[role=dialog]")')


def archive(b):
    path=DIRECTORY/'collagene-workspace.cg2workspace';path.unlink(missing_ok=True)
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.click('Workspace');b.click('Save workspace file');value=wait_file(path);close(b)
    return manifest(value)


def oracle(report,sequence):
    """Independent original-coordinate/IUPAC oracle; no application code imported."""
    profiles={p['enzyme']:p for p in report['profiles']}
    n=len(sequence)
    assert report['source_length']==n
    for site in report['sites']:
        pattern=profiles[site['enzyme']]['recognition']
        if site['reverse']:pattern=reverse(pattern)
        positions=[(site['recognition_start']+i)%n for i in range(len(pattern))]
        observed=''.join(sequence[i] for i in positions)
        assert all(MASK[a]&MASK[b] for a,b in zip(observed,pattern))
        assert site['certain']==all(a in 'ACGT' for a in observed)
        assert len(site['effects'])==3
        for effect in site['effects']:
            for context in effect['overlaps']:
                expected=[(context['start']+i)%n for i in range(len(context['motif']))]
                assert context['positions']==expected
                seq=''.join(sequence[i] for i in expected)
                assert all(MASK[a]&MASK[b] for a,b in zip(seq,context['motif']))
                assert context['certain']==all(a in 'ACGT' for a in seq)
                for base in context['modified_bases']:
                    assert base['position'] in expected
                    actual=sequence[base['position']]
                    if base['reverse']:actual=reverse(actual)
                    assert MASK[actual]&MASK[base['base']],(actual,base)


def single(r,enzyme):
    values=[s for s in r['sites'] if s['enzyme']==enzyme]
    assert len(values)==1,values
    return values[0]


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for p in DIRECTORY.iterdir():
        if p.is_file():p.unlink()
    b=Browser();b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':ORIGIN});b.wait('document.body?.innerText.includes("Saved in this browser")')
    fixtures=[('dam_sensitive','AAATCTAGATCAAA',False),('cpg_sensitive','AAAGCGGCCGCAAACCCGGGAAA',False),('some_overlap','CGATATCG',False),('unsupported_formulation','CCAGGTCTCAAAACGTCTCAAA',False),('ambiguous_recognition','AAATCTAGNCAAA',False),('reverse_type_iis','AAAGAGACCAAA',False),('circular_boundary','GATCAAATCTA',True),('many_sites','TCTAGATCAAA'*180,False),('dam_sensitive','AAATCTAGATCAAA',False)]
    paste(b,''.join(record(*fixture)for fixture in fixtures))
    active(b,'dam_sensitive');panel(b)
    assert b.evaluate('Array.from(document.querySelectorAll("select[aria-label$=assumption]")).every(e=>e.value==="unknown")')
    b.evaluate('document.querySelector(".methylation-enzyme-picker").open=true')
    assert b.evaluate('document.querySelectorAll(".methylation-enzyme-choice").length')==12
    assert b.evaluate('(()=>{const row=document.querySelector(".methylation-enzyme-choice"),input=row.querySelector("input"),label=row.querySelector("span");return getComputedStyle(row).flexDirection==="row"&&input.getBoundingClientRect().width<32&&label.getBoundingClientRect().left>=input.getBoundingClientRect().right})()'),'Profile selection rows must retain horizontal checkbox layout'
    select(b,['XbaI']);states(b,dam='unknown');assess(b)
    unknown=report(b);r=unknown['report'];assert r['complete'];oracle(r,fixtures[0][1]);assert single(r,'XbaI')['outcome']=='Unknown'
    before=archive(b)
    states(b,dam='fully');assess(b);full=report(b);r=full['report'];site=single(r,'XbaI');assert site['outcome']=='Blocked';oracle(r,fixtures[0][1])
    dam=next(e for e in site['effects']if e['category']=='Dam');assert dam['overlaps'] and any(o['certain']for o in dam['overlaps'])
    # Exact site/context positions and modeled target strands are inspectable.
    b.evaluate('document.querySelector(".methylation-site-entry").click()');b.wait('!!document.querySelector(".methylation-site-detail")')
    assert 'original source bases' in b.evaluate('document.querySelector(".methylation-site-detail").textContent')
    detail=b.evaluate('document.querySelector(".methylation-site-detail").textContent')
    assert 'Some(' not in detail and 'BlockedByOverlapping' not in detail
    assert f"Top cut boundary: {site['top_cut']} bases from source start" in detail
    assert 'Blocked by an overlapping motif' in detail
    b.click('Open target base');b.wait('!document.querySelector(".methylation-panel")')
    target=dam['overlaps'][0]['modified_bases'][0]['position']+1
    b.wait('document.querySelector(".statusbar")?.textContent.includes('+json.dumps(f'Selection {target}–{target}')+')')
    # The exact source event DAG is unchanged: review is not a molecule edit.
    assert archive(b)['documents']==before['documents']
    panel(b);select(b,['XbaI']);states(b);assess(b);assert single(report(b)['report'],'XbaI')['outcome']=='PermittedUnderProfile'
    states(b,dam='mixed');assess(b);assert single(report(b)['report'],'XbaI')['outcome']=='Unknown'
    # CpG sensitivity uses actual source motifs, not methylation guessed from a name.
    active(b,'cpg_sensitive');panel(b);select(b,['NotI','SmaI']);states(b,cpg='fully');assess(b);cp=report(b)['report'];oracle(cp,fixtures[1][1])
    assert single(cp,'NotI')['outcome']=='Blocked' and single(cp,'SmaI')['outcome']=='Blocked'
    # Some-combination rules and formulation conflicts stay unresolved.
    active(b,'some_overlap');panel(b);select(b,['EcoRV']);states(b,cpg='fully');assess(b);sr=report(b)['report'];oracle(sr,fixtures[2][1]);assert single(sr,'EcoRV')['outcome']=='Unknown'
    active(b,'unsupported_formulation');panel(b);select(b,['BsaI','BsmBI']);states(b,dam='fully',dcm='fully',cpg='fully');assess(b);ur=report(b)['report'];oracle(ur,fixtures[3][1])
    assert single(ur,'BsaI')['outcome']=='Unknown' and single(ur,'BsmBI')['outcome']=='Unknown'
    assert any('R0535' in p['supplier_identity']for p in ur['profiles'])
    assert any('not silently' in note.lower() for p in ur['profiles']for note in p['notes'])
    active(b,'ambiguous_recognition');panel(b);select(b,['XbaI']);states(b);assess(b);ar=report(b)['report'];oracle(ar,fixtures[4][1]);assert any(not s['certain'] and s['outcome']=='Unknown'for s in ar['sites'])
    active(b,'reverse_type_iis');panel(b);select(b,['BsaI']);states(b);assess(b);rr=report(b)['report'];oracle(rr,fixtures[5][1]);assert single(rr,'BsaI')['reverse'];assert single(rr,'BsaI')['top_cut'] is None or single(rr,'BsaI')['bottom_cut'] is None
    b.evaluate('document.querySelector(".methylation-site-entry").click()')
    b.wait('!!document.querySelector(".methylation-site-detail")')
    assert 'unavailable in this source' in b.evaluate('document.querySelector(".methylation-site-detail").textContent')
    active(b,'circular_boundary');panel(b);select(b,['XbaI']);states(b,dam='fully');assess(b);cr=report(b)['report'];oracle(cr,fixtures[6][1]);cs=single(cr,'XbaI');assert cs['recognition_start']+cs['recognition_length']>cr['source_length'];assert cs['outcome']=='Blocked'
    # Resource exhaustion has explicit incomplete evidence and cannot establish absent sites.
    active(b,'many_sites');panel(b);select(b,['XbaI']);states(b,dam='fully');fill(b,'Methylation site limit','1');assess(b);lr=report(b)['report'];assert not lr['complete'] and lr['limits'] and len(lr['sites'])==1
    states(b,dam='fully');fill(b,'Methylation comparison budget','1');assess(b);lr=report(b)['report'];assert not lr['complete'] and lr['limits']
    states(b,dam='fully')
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Assess methylation").click();queueMicrotask(()=>Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Cancel methylation assessment").click());})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".methylation-result")===null')
    assert 'cancelled' in b.evaluate('document.querySelector(".methylation-status").textContent').lower()
    states(b,dam='fully');assess(b);large=report(b)['report'];assert large['complete'] and len(large['sites'])==180
    assert b.evaluate('document.querySelectorAll(".methylation-site-entry").length')==25
    b.click('Next methylation sites');assert b.evaluate('Number(document.querySelector(".methylation-site-entry").dataset.siteIndex)')==25
    fill(b,'Methylation site outcome','Blocked','change');assert b.evaluate('document.querySelectorAll(".methylation-site-entry").length')==25
    b.click('Expand methylation');b.wait('document.querySelector(".methylation-panel").classList.contains("methylation-expanded")')
    b.evaluate('document.querySelector(".methylation-site-entry").click()')
    b.wait('!!document.querySelector(".methylation-site-detail")')
    b.evaluate('document.querySelector(".methylation-site-detail").scrollIntoView({block:"center"})')
    b.screenshot('/tmp/collagene2-methylation-desktop.png')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        b.evaluate('document.querySelector(".methylation-site-detail").scrollIntoView({block:"center"})')
        assert b.evaluate('document.documentElement.scrollWidth<=innerWidth')
        b.screenshot('/tmp/collagene2-methylation-mobile.png')
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False});b.click('Close expanded methylation')
    # Option and identity changes cancel pending work before it can become current.
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Assess methylation").click();const e=document.querySelector("[aria-label=\\"Methylation Dam assumption\\"]");e.value="unmodified";e.dispatchEvent(new Event("change",{bubbles:true}));})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".methylation-result")===null')
    active(b,'dam_sensitive');panel(b);select(b,['XbaI']);states(b,dam='fully');assess(b);old=report(b)
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Assess methylation").click();Array.from(document.querySelectorAll(".document-entry")).filter(e=>e.querySelector("strong").textContent==="dam_sensitive").at(-1).click();})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".methylation-result")===null')
    assess(b);new=report(b);assert new['source_document_id']!=old['source_document_id'] and new['source_revision']==old['source_revision']
    # Unmount must stop the pending job, and returns with no trusted cache.
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Assess methylation").click();Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Features").click();})()')
    time.sleep(.2);panel(b);assert b.evaluate('document.querySelector(".methylation-result")===null')
    # The virtual digest retains cleavage assumptions and its visible disclosure
    # bridge back to this separate planning review.
    select(b,['XbaI']);states(b,dam='fully');assess(b)
    b.click('Review cleavage assumptions');b.wait('!document.querySelector(".methylation-panel")')
    assert 'Products assume complete cleavage; methylation, reaction conditions and laboratory success are not assessed.' in b.evaluate('document.querySelector(".digest-panel").innerText')
    b.pointer_click('document.querySelector(".digest-options > summary")')
    b.wait('document.querySelector(".digest-options").open')
    assert 'Review methylation assumptions' in b.evaluate('document.querySelector(".digest-options").innerText')
    b.pointer_click('Array.from(document.querySelectorAll(".digest-options button")).find(e=>e.textContent==="Review methylation assumptions")')
    b.wait('!!document.querySelector(".methylation-panel")')
    assert b.evaluate('document.querySelector(".methylation-result")===null')
    # No saved report, physical molecule metadata, products or event changes.
    final=archive(b);assert final['documents']==before['documents'] and final['history']==before['history']
    assert not any(a['kind'].startswith('methylation')for a in final['attachments'].values())
    paste(b,'>methylation_non_dna [molecule=protein]\nMPEPTIDE')
    panel(b)
    assert b.evaluate('Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Assess methylation").disabled')
    assert b.evaluate('document.querySelector(".methylation-result")===null')
    assert not b.exceptions and not b.console_errors,(b.exceptions,b.console_errors)
    print('PASS: methylation review explicit states, XbaI Dam/CpG blocks, unresolved profile/formulation cases, independent circular/reverse/IUPAC coordinates, exact identity/options cancellation, original-target navigation, bounded evidence/mobile and unchanged source history'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__=='__main__':main()
