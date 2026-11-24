"""Literal TA/GC chemistry and oriented overlap candidates in the actual WASM UI.
Uses only synthetic sequences; external pdftotext is an independent dev validator.
"""
import json
import subprocess
import time
from pathlib import Path
from browser_qa import QA_ORIGIN, Browser, mobile_enabled
from workspace_browser_qa import wait_file,manifest

DIRECTORY=Path('/tmp/collagene2-methods-qa')


def record(name,sequence,circular=False):
    return f'LOCUS       {name} {len(sequence)} bp DNA '+('circular'if circular else'linear')+'\nFEATURES             Location/Qualifiers\nORIGIN\n        1 '+sequence.lower()+'\n//\n'


def fill(b,label,value,event='input'):
    selector='[aria-label='+json.dumps(label)+']'
    b.evaluate('(()=>{const e=document.querySelector('+json.dumps(selector)+');e.value='+json.dumps(value)+';e.dispatchEvent(new Event('+json.dumps(event)+',{bubbles:true}));})()')


def check(b,label,enabled=True):
    selector='input[aria-label='+json.dumps(label)+']'
    b.evaluate('(()=>{const e=document.querySelector('+json.dumps(selector)+');if(e.checked!=='+json.dumps(enabled)+')e.click();})()')


def source(b,label,name):
    selector='select[aria-label='+json.dumps(label)+']'
    b.evaluate('(()=>{const e=document.querySelector('+json.dumps(selector)+');const o=Array.from(e.options).find(o=>o.textContent==='+json.dumps(name)+');if(!o)throw Error("Missing source");e.value=o.value;e.dispatchEvent(new Event("change",{bubbles:true}));})()')


def active(b,name):
    b.evaluate('Array.from(document.querySelectorAll(".document-entry")).find(e=>e.querySelector("strong")?.textContent==='+json.dumps(name)+').click()')
    b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name))


def methods(b):
    b.click('Cloning methods');b.wait('!!document.querySelector(".cloning-methods-panel")')


def compute(b):
    b.click('Compute method evidence');b.wait('!!document.querySelector(".method-evidence")')


def verify(b,status):
    b.click('Review method candidate');b.wait('!!document.querySelector(".candidate-review")')
    b.click('Verify candidate');b.wait('document.querySelector(".verification-result strong")?.textContent==='+json.dumps(status))


def evidence(b):
    p=DIRECTORY/'cloning-method-evidence.json';p.unlink(missing_ok=True)
    b.click('Download method evidence');return json.loads(wait_file(p))


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for p in DIRECTORY.iterdir():
        if p.is_file():p.unlink()
    b=Browser();b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':QA_ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':QA_ORIGIN});b.wait('document.body?.innerText.includes("Saved in this browser")')
    # EcoRV cuts AT[GAT|ATC]GG at reference5. Rotated blunt vector is
    # ATCGGATGAT. Two ordered T/A junctions retain BOTH added reference bases.
    fixtures=[('ta_vector','ATGATATCGG',True),('pcr_template','ACGACCCCGGGGTAAC',False),('ta_target','ATCGGATGATTACGACCCCGGGGTAACA',True),('gc_target','ATCGGATGATCACGACCCCGGGGTAACG',True),('missing_tails','ATCGGATGATACGACCCCGGGGTAAC',True),('ov_left','AACCTTGCAAACGTAGCTACGATTC',False),('ov_right_rc','TAAGGATCCGAATCGTAGCTACGT',False),('ov_target','AACCTTGCAAACGTAGCTACGATTCGGATCCTTA',False),('short_left','AAAACCCC',False),('short_right','CCCCGGGG',False),('short_target','AAAACCCCGGGG',False)]
    b.click('Paste sequence');fill(b,'Sequence input',''.join(record(*f)for f in fixtures));b.click('Open sequence')
    b.wait('document.querySelector("h1")?.textContent==="ta_vector"')
    # The final document remains pointer-accessible inside a genuinely scrollable
    # library; programmatic .click() alone would hide a clipped-sidebar bug.
    b.evaluate('document.querySelector(".document-entry:last-of-type")?.scrollIntoView({block:"nearest"})')
    position=b.evaluate('(()=>{const list=Array.from(document.querySelectorAll(".document-entry")),button=list[list.length-1],r=button.getBoundingClientRect(),library=document.querySelector(".library");button.scrollIntoView({block:"nearest"});const q=button.getBoundingClientRect();return {x:q.left+q.width/2,y:q.top+q.height/2,top:q.top,bottom:q.bottom,scroll:library.scrollTop,overflow:library.scrollHeight>library.clientHeight};})()')
    assert position['overflow'] and position['scroll']>0 and 72<=position['top']<position['bottom']<=1000,position
    b.call('Input.dispatchMouseEvent',{'type':'mousePressed','x':position['x'],'y':position['y'],'button':'left','clickCount':1})
    b.call('Input.dispatchMouseEvent',{'type':'mouseReleased','x':position['x'],'y':position['y'],'button':'left','clickCount':1})
    b.wait('document.querySelector("h1")?.textContent==="short_target"')
    active(b,'ta_target');methods(b)
    source(b,'Vector tailing source','ta_vector');source(b,'Insert tailing source','pcr_template')
    fill(b,'Insert tailing forward primer','ACGA');fill(b,'Insert tailing reverse primer','GTTA')
    compute(b);r=evidence(b);assert len(r)==2
    assert r[0]['tailing']['source']['top']['bases']=='ATCGGATGAT'
    assert r[0]['tailing']['duplex']['top']['bases']=='ATCGGATGATT'
    assert r[0]['tailing']['duplex']['bottom']['bases']=='ATCATCCGATT'
    assert r[0]['tailing']['duplex']['bottom_start']==-1
    assert r[1]['tailing']['duplex']['top']['bases']=='ACGACCCCGGGGTAACA'
    assert r[1]['tailing']['duplex']['left']['five_prime']=='Hydroxyl'
    assert all(len(p['tailing']['additions'])==2 for p in r)
    assert all(not m['synthesized'] for p in r for m in p['tailing']['mappings'])
    # PCR5′ OH cannot close the two insert-donor bonds without another step.
    verify(b,'Invalid')
    assert 'phosphate' in b.evaluate('document.querySelector(".verification-result").textContent').lower()
    assert 'Reconstructed product matches' not in b.evaluate('document.querySelector(".verification-result").textContent')
    # Configuration survives candidate review; phosphorylation is explicit graphdata.
    methods(b);assert b.evaluate('document.querySelector('+json.dumps('[aria-label="Insert tailing forward primer"]')+').value')=='ACGA'
    fill(b,'Insert tailing terminal treatment','phosphorylate','change');compute(b);verify(b,'ConsistentUnderAssumptions')
    text=b.evaluate('document.querySelector(".verification-result").textContent')
    assert '28 nt' in text and 'Covalently closed in this model: true' in text,text
    b.click('Save candidate');assert wait_file(DIRECTORY/'candidate.cg2plan')[4:6]==bytes([7,0])
    b.click('Export HTML protocol');html=wait_file(DIRECTORY/'cloning-protocol.html').decode();assert 'ACGA' in html and 'GTTA' in html and 'tail' in html.lower()
    b.click('Export PDF protocol');pdf=DIRECTORY/'cloning-protocol.pdf';wait_file(pdf)
    extracted=subprocess.check_output(['pdftotext',str(pdf),'-'],text=True)
    assert 'ACGA' in extracted and 'GTTA' in extracted and 'tail' in extracted.lower()
    # The candidate stays bound to its independent target. Missing the two tail
    # bases is a different target, not a reason to reinterpret the plan asverified.
    active(b,'missing_tails');b.click('Verify candidate');b.wait('document.querySelector(".notice")?.textContent.toLowerCase().includes("target")');assert not b.evaluate('document.querySelector(".verification-result")')
    active(b,'gc_target');methods(b);b.click('Use C/G pair')
    fill(b,'Vector tailing terminal treatment','dephosphorylate','change');compute(b);verify(b,'Invalid')
    methods(b);check(b,'Allow retained tailing ligation nicks');compute(b);verify(b,'ConsistentUnderAssumptions')
    text=b.evaluate('document.querySelector(".verification-result").textContent')
    assert 'Covalently closed in this model: false' in text and 'Retained nicks' in text,text
    # Nativearchives retain the new operation semantics underplan_schema3.
    b.click('Workspace');b.click('Save workspace file');archive=wait_file(DIRECTORY/'collagene-workspace.cg2workspace');assert manifest(archive)['plan_schema']==5
    b.evaluate('document.querySelector('+json.dumps('button[aria-label="Close dialog"]')+').click()')
    # Reverse-complement the supplied right source explicitly. The checker gives
    # ORIGINAL-reference coordinates, then verifier independently joins theDNA.
    active(b,'ov_target');methods(b);fill(b,'Cloning workflow','overlap','change');check(b,'Method product circular',False)
    source(b,'Overlap source 1','ov_left');source(b,'Overlap source 2','ov_right_rc');check(b,'Reverse overlap source 2')
    fill(b,'Profile overlap lengths','15');compute(b);r=evidence(b)
    assert r['geometry_valid'] and r['complete'] and r['recommendation_compliant']
    assert r['orientations']==[False,True]
    assert r['junctions'][0]['homology']=='ACGTAGCTACGATTC'
    assert r['junctions'][0]['downstream_homology']=='ACGTAGCTACGATTC'
    assert r['junctions'][0]['upstream_span']=={'start':10,'end':25,'reverse':False}
    assert r['junctions'][0]['downstream_span']=={'start':9,'end':24,'reverse':True}
    assert not r['exhaustive_alternative_search']
    verify(b,'ConsistentUnderAssumptions')
    # A bounded incomplete repeat search cannot borrow an earlier full report.
    methods(b);fill(b,'Overlap repeat budget','20');b.wait('!document.querySelector(".method-evidence")')
    compute(b);r=evidence(b);assert r['geometry_valid'] and not r['complete'];verify(b,'Incomplete')
    # Recommended lengths are advisory. Exact4bp sequencegeometry remainsdistinct
    # from whether the reaction meets the documented HiFi guidance.
    active(b,'short_target');methods(b);fill(b,'Overlap repeat budget','5000000')
    source(b,'Overlap source 1','short_left');source(b,'Overlap source 2','short_right');check(b,'Reverse overlap source 2',False);fill(b,'Profile overlap lengths','4')
    compute(b);r=evidence(b);assert r['geometry_valid'] and r['complete'] and r['recommendation_compliant'] is False
    assert 'outside published guidance' in b.evaluate('document.querySelector(".overlap-profile-evidence").textContent')
    verify(b,'ConsistentUnderAssumptions')
    # Incompatible geometry cannot become a valid target just by selecting a brand.
    methods(b);fill(b,'Profile overlap lengths','5');fill(b,'Overlap profile method','InFusion','change');compute(b);r=evidence(b);assert not r['geometry_valid'];verify(b,'Invalid')
    # Source/options changes cancel a pending job before it can restore oldevidence.
    methods(b);fill(b,'Profile overlap lengths','4')
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(b=>b.textContent.trim()==="Compute method evidence").click();const e=document.querySelector('+json.dumps('[aria-label="Profile overlap lengths"]')+');e.value="2";e.dispatchEvent(new Event("input",{bubbles:true}));})()')
    time.sleep(.5);assert not b.evaluate('document.querySelector(".method-evidence")')
    compute(b);assert not b.evaluate('document.querySelector(".cloning-methods-panel").textContent.includes("on : click")')
    b.screenshot('/tmp/collagene2-cloning-methods.png')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        b.wait('document.querySelector(".cloning-methods-panel")')
        assert b.evaluate('document.documentElement.scrollWidth<=390')
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    print('PASS: literal TA/GCjunctionbases, PCRphosphates/sealed-vs-nickedproducts, independenttarget rejection, newcandidate/workspaceversions, actualHTML/PDFoligos, orientedsourcehomologycoordinates, advisoryprofiles, incompletebudgets, stalejobcancellation andmobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))

if __name__=='__main__':main()
