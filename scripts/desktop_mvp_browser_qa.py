#!/usr/bin/env python3
"""Desktop viewer/GenBank/virtual-digest acceptance with independent literal oracles."""
import base64
import hashlib
import io
import json
import random
import re
import time
import zipfile
from pathlib import Path
from browser_qa import Browser
from methods_browser_qa import fill
from workspace_browser_qa import wait_file,slots

DIRECTORY=Path('/tmp/collagene2-desktop-mvp-qa')
ORIGIN='http://127.0.0.1:8080'
ENZYMES={'EcoRI':('GAATTC',1,5),'BamHI':('GGATCC',1,5),'NotI':('GCGGCCGC',2,6),'SapI':('GCTCTTC',8,11)}
COMP=str.maketrans('ACGT','TGCA')


def rc(s):return s.translate(COMP)[::-1]


def record(name,sequence,circular=False,features=''):
    lines='\n'.join(f'{i+1:9d} '+sequence[i:i+60].lower()for i in range(0,len(sequence),60))
    return f'LOCUS       {name:<16} {len(sequence):>10} bp    DNA     {"circular"if circular else "linear  "} SYN 01-JAN-2026\nDEFINITION  Independently authored desktop acceptance fixture.\nACCESSION   AB123456\nVERSION     AB123456.1\nCOMMENT     Original fixture evidence: DESKTOP_MVP_LITERAL.\nFEATURES             Location/Qualifiers\n{features}ORIGIN\n{lines}\n//\n'


def upload(b,path,selector='input[type=file]'):
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId']
    node=b.call('DOM.querySelector',{'nodeId':root,'selector':selector})['nodeId'];assert node,selector
    b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})


def click_expression(b,expression):
    rect=b.evaluate('(()=>{const e='+expression+';if(!e)throw Error("Missing acceptance target");e.scrollIntoView({block:"center",inline:"nearest"});const r=e.getBoundingClientRect();return {x:r.x+r.width/2,y:r.y+r.height/2}})()')
    assert b.evaluate('(()=>{const e='+expression+',h=document.elementFromPoint('+str(rect['x'])+','+str(rect['y'])+');return e===h||e.contains(h)})()'),'Target is obscured'
    b.call('Input.dispatchMouseEvent',{'type':'mouseMoved',**rect})
    b.call('Input.dispatchMouseEvent',{'type':'mousePressed','button':'left','clickCount':1,**rect})
    b.call('Input.dispatchMouseEvent',{'type':'mouseReleased','button':'left','clickCount':1,**rect})


def document(b,name,index=0):
    expression='Array.from(document.querySelectorAll(".document-entry")).filter(e=>e.querySelector("strong")?.textContent==='+json.dumps(name)+')['+str(index)+']'
    click_expression(b,expression);b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name))


def download(b,label,filename):
    path=DIRECTORY/filename;path.unlink(missing_ok=True);click(b,label);return wait_file(path)


def click(b,label):
    text=json.dumps(label)
    if b.evaluate('Array.from(document.querySelectorAll("[data-tool-id],[data-panel-id]")).some(e=>e.textContent.trim()==='+text+')'):
        b.click(label);return
    target='Array.from(document.querySelectorAll("button")).find(e=>e.textContent.trim()==='+text+')'
    assert b.evaluate('!!('+target+')'),label
    ancestor='(()=>{const closed=[];let p=('+target+').parentElement;while(p){if(p.tagName==="DETAILS"&&!p.open)closed.push(p);p=p.parentElement;}return closed.at(-1)})()'
    while b.evaluate('!!('+ancestor+')'):
        b.pointer_click('('+ancestor+').querySelector("summary")')
    b.pointer_click(target)


def set_enzymes(b,names):
    for name in ['EcoRI','BamHI','HindIII','EcoRV','SmaI','NotI','XbaI','BsaI','PstI','BsmBI','BbsI','SapI']:
        selector='input[aria-label='+json.dumps('Digest enzyme '+name)+']'
        target='document.querySelector('+json.dumps(selector)+')'
        if b.evaluate('('+target+').checked')!=(name in names):b.pointer_click(target)


def digest(b,names):
    click(b,'Digest');b.wait('!!document.querySelector(".digest-panel")');set_enzymes(b,names)
    click(b,'Run virtual digest');b.wait('!!document.querySelector(".digest-results")',timeout=60)
    return json.loads(download(b,'Download digest evidence','virtual-digest.json'))


def oracle(sequence,circular,names):
    n=len(sequence);sites=[];cuts=set()
    for name in names:
        motif,top,bottom=ENZYMES[name]
        for rev in [False,True]:
            pattern=rc(motif)if rev else motif
            if rev and pattern==motif:continue
            for start in range(n if circular else max(0,n-len(motif)+1)):
                text=''.join(sequence[(start+i)%n]for i in range(len(motif)))
                if text!=pattern:continue
                t=start+(len(motif)-bottom if rev else top);u=start+(len(motif)-top if rev else bottom)
                valid=circular or(0<=t<=n and 0<=u<=n)
                sites.append((name,start,rev,t%n if circular else t,u%n if circular else u,valid))
                if valid:cuts.add(t%n if circular else t)
    if not all(s[-1]for s in sites):return sites,None
    if not cuts:return sites,[n]
    positions=sorted(cuts);lengths=[]
    if not circular:positions=[0]+positions+[n]
    lengths=[z-a for a,z in zip(positions,positions[1:])if z>a]
    if circular:lengths.append(n-positions[-1]+positions[0])
    return sites,sorted(lengths,reverse=True)


def gb_sequence(text):
    body=text.split('ORIGIN',1)[1].split('//',1)[0]
    return ''.join(c for c in body if c.isalpha()).upper()


def digest_comment(text,prefix):
    chunks=[];total=None
    for line in text.splitlines():
        match=re.search(re.escape(prefix)+r' (\d+)/(\d+) ([A-Za-z0-9+/=]+)$',line)
        if not match:continue
        index,count,encoded=match.groups();assert int(index)==len(chunks)+1
        if total is None:total=int(count)
        assert int(count)==total;chunks.append(encoded)
    assert chunks and len(chunks)==total
    return json.loads(base64.b64decode(''.join(chunks),validate=True))


def capture_map(b):
    rect=b.evaluate('(()=>{const r=document.querySelector(".map-stage").getBoundingClientRect();return {x:r.x,y:r.y,width:r.width,height:r.height,scale:1}})()')
    return hashlib.sha256(base64.b64decode(b.call('Page.captureScreenshot',{'format':'png','clip':rect})['data'])).hexdigest()


def hold_jobs(b):
    b.evaluate('window.__heldPosts=[];window.__holdJobs=true;window.__post=Worker.prototype.postMessage;Worker.prototype.postMessage=function(...args){if(window.__holdJobs){window.__heldPosts.push([this,args]);return;}return window.__post.apply(this,args)}')


def release_jobs(b):
    b.evaluate('(()=>{window.__holdJobs=false;const pending=window.__heldPosts.splice(0);Worker.prototype.postMessage=window.__post;for(const [worker,args]of pending){try{window.__post.apply(worker,args)}catch(error){}}})()')


def archive(b):
    b.wait('document.body?.innerText.includes("Saved in this browser")');click(b,'Workspace')
    data=download(b,'Save workspace file','collagene-workspace.cg2workspace')
    b.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()');return data


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for path in DIRECTORY.iterdir():
        if path.is_file():path.unlink()
    random.seed(23);dna=list(''.join(random.choice('ACGT')for _ in range(120)))
    for start,motif in [(14,'GAATTC'),(45,'GGATCC'),(70,'GAATTC'),(100,'GCGGCCGC')]:dna[start:start+len(motif)]=motif
    dna[:3]='TTC';dna[-3:]='GAA';dna=''.join(dna)
    features='     source          1..120\n                     /organism="synthetic construct"\n     CDS             complement(join(111..120,1..12))\n                     /label="Origin reverse feature"\n     promoter        12..24\n                     /label="Promoter marker"\n     misc_feature    26..44\n                     /label="Retained middle"\n'
    fixtures=[('MVP_CIRCULAR',dna,True,features),('MVP_LINEAR',dna,False,'     misc_feature    26..44\n                     /label="Linear marker"\n'),('MVP_EDGE','GCTCTTCA',False,''),('MVP_AMBIG','AAAAGAANTCAAAA',False,''),('MVP_UNCUT','A'*47,True,''),('MVP_MANY',('GAATTC'+'ACGTAC')*60,True,''),('pCG_demo','ACGTACGT',False,'')]
    paths={}
    for name,seq,circ,feature in fixtures:
        path=DIRECTORY/(name+'.gb');path.write_text(record(name,seq,circ,feature));paths[name]=path
    b=Browser();b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':ORIGIN,'storageTypes':'indexeddb'});b.exceptions.clear();b.console_errors.clear()
    b.call('Page.navigate',{'url':ORIGIN});b.wait('document.body?.innerText.includes("Saved in this browser")')
    assert b.evaluate('Array.from(document.querySelectorAll(".mvp-inspector-tabs > button")).map(e=>e.textContent).join(",")==="Features,Digest"')
    assert b.evaluate('!!Array.from(document.querySelectorAll(".top-actions button")).find(e=>e.textContent==="Open GenBank")')
    assert b.evaluate('!!document.querySelector(".mvp-example")')
    # Error/size checks are atomic and do not replace the currently open example.
    bad=DIRECTORY/'bad.gb';bad.write_text('LOCUS broken\nnot a GenBank record\n');before=b.evaluate('document.querySelectorAll(".document-entry").length');upload(b,bad);b.wait('!!document.querySelector(".notice")');assert b.evaluate('document.querySelectorAll(".document-entry").length')==before
    oversized=DIRECTORY/'oversized.gb'
    with oversized.open('wb')as f:f.truncate(33*1024*1024)
    upload(b,oversized);b.wait('document.querySelector(".notice")?.textContent.includes("32 MiB")');assert b.evaluate('document.querySelectorAll(".document-entry").length')==before
    # Import actual annotated GenBank; the examples remain clearly separate.
    upload(b,paths['MVP_CIRCULAR']);b.wait('document.querySelector("h1")?.textContent==="MVP_CIRCULAR"')
    assert not b.evaluate('!!document.querySelector(".mvp-example")')
    assert b.evaluate('document.querySelector(".map-panel h2").textContent==="Plasmid map"')
    assert b.evaluate('!!document.querySelector(".map-stage canvas").getContext("webgl2")')
    source_native=download(b,'Save .cg2','MVP_CIRCULAR.cg2');assert paths['MVP_CIRCULAR'].read_bytes()in source_native
    old_image=capture_map(b);rect=b.evaluate('(()=>{const r=document.querySelector(".map-stage canvas").getBoundingClientRect();return {x:r.x+r.width/2,y:r.y+r.height/2,r:Math.min(r.width,r.height)*.34}})()')
    b.call('Input.dispatchMouseEvent',{'type':'mousePressed','button':'left','buttons':1,'clickCount':1,'x':rect['x']+rect['r'],'y':rect['y']})
    for dx,dy in [(0.9,.4),(.7,.7),(.4,.9),(0,1)]:b.call('Input.dispatchMouseEvent',{'type':'mouseMoved','button':'left','buttons':1,'x':rect['x']+rect['r']*dx,'y':rect['y']+rect['r']*dy})
    b.call('Input.dispatchMouseEvent',{'type':'mouseReleased','button':'left','clickCount':1,'x':rect['x'],'y':rect['y']+rect['r']})
    b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))');assert capture_map(b)!=old_image
    assert download(b,'Save .cg2','MVP_CIRCULAR.cg2')==source_native,'Rotating changed scientific source data'
    click(b,'Features');click_expression(b,'Array.from(document.querySelectorAll(".feature-entry")).find(e=>e.querySelector("strong").textContent==="Origin reverse feature")')
    b.wait('document.querySelector(".annotation-sequence")?.textContent==='+json.dumps(rc(dna[110:]+dna[:12])))
    assert 'complement(join(111..120,1..12))'in b.evaluate('document.querySelector(".annotation-details").textContent')
    # Contiguous sequence selection uses actual mouse events and preserves source data.
    region=b.evaluate('(()=>{const e=document.querySelector(".sequence-scroll");e.scrollTop=0;e.dispatchEvent(new Event("scroll"));const r=document.querySelector(".sequence-canvas").getBoundingClientRect();return {x:r.x,y:r.y}})()')
    for col,shift in [(0,False),(7,True)]:
        point={'x':region['x']+68+col*8.5,'y':region['y']+25,'modifiers':8 if shift else 0}
        b.call('Input.dispatchMouseEvent',{'type':'mousePressed','button':'left','clickCount':1,**point});b.call('Input.dispatchMouseEvent',{'type':'mouseReleased','button':'left','clickCount':1,**point})
    b.wait('document.querySelector(".statusbar")?.textContent.includes("1–8")')
    evidence=digest(b,['EcoRI','BamHI']);expected_sites,expected_lengths=oracle(dna,True,['EcoRI','BamHI']);report=evidence['report']
    assert report['scan_complete'] and report['products_determined'] and report['ref_top_lengths']==expected_lengths
    actual={(s['enzyme'],s['recognition_start'],s['reverse'],s['top_cut'],s['bottom_cut'])for s in report['sites']}
    assert actual=={s[:5]for s in expected_sites}
    wrapping=next(s for s in report['sites']if len(s['spans'])==2);assert sum(s['end']-s['start']for s in wrapping['spans'])==6
    click_expression(b,'document.querySelectorAll(".digest-site button")['+str(wrapping['index'])+']')
    assert b.evaluate('document.querySelector(".statusbar").textContent.includes("Click a base")'),'Origin-spanning sites must not create a misleading contiguous edit range'
    assert report['fragments'][0]['top_start']==15 and report['fragments'][0]['duplex']['top_length']==31
    ends=report['fragments'][0]['duplex']['ends'];assert ends[0]['polarity']=='FivePrime'and ends[0]['overhang']=='AATT'and ends[1]['polarity']=='FivePrime'and ends[1]['overhang']=='GATC'
    product=download(b,'Export fragment GenBank','digest-fragment-1.gb').decode();assert gb_sequence(product)==dna[15:46]
    assert 'Retained middle'in product and '11..29'in product and 'collagene_review_required'in product
    assert 'unverified'in product.lower()
    saved_recipe=digest_comment(product,'CG-DIGEST-RECIPE-v1');saved_ends=digest_comment(product,'CG-DIGEST-ENDS-v1')
    assert saved_recipe['inspection_revision']==report['batch_revision'] and saved_recipe['source_role']==report['source_role'] and saved_recipe['fragment_index']==0
    assert saved_ends['top_length']==31 and saved_ends['bottom_length']==31 and saved_ends['ends']==report['fragments'][0]['duplex']['ends']
    # Native installation is atomic, including a synchronous second-put failure.
    b.wait('document.body?.innerText.includes("Saved in this browser")');before_slots=slots(b);before_documents=b.evaluate('document.querySelectorAll(".document-entry").length')
    b.evaluate('window.__put=IDBObjectStore.prototype.put;IDBObjectStore.prototype.put=function(value,key){if(this.transaction.db.name==="collagene"&&this.name==="workspace"&&key==="current")throw new DOMException("Synthetic second-put quota failure","QuotaExceededError");return window.__put.call(this,value,key)}')
    click(b,'Open fragment');b.wait('document.querySelector(".digest-open-panel")?.textContent.includes("fragment was not opened")',timeout=60)
    b.evaluate('IDBObjectStore.prototype.put=window.__put');assert slots(b)==before_slots and b.evaluate('document.querySelectorAll(".document-entry").length')==before_documents and b.evaluate('document.querySelector("h1").textContent')=='MVP_CIRCULAR'
    b.evaluate('document.querySelector('+json.dumps('[aria-label="Close dialog"]')+').click()');b.wait('!document.querySelector("[role=dialog]")')
    # Close the modal in a microtask after the real current put is queued, before
    # IDB can commit. No production-only hook or synthetic success receipt exists.
    before_cancel=slots(b)
    b.evaluate('window.__put=IDBObjectStore.prototype.put;window.__cancelledPending=false;IDBObjectStore.prototype.put=function(value,key){const request=window.__put.call(this,value,key);if(this.transaction.db.name==="collagene"&&this.name==="workspace"&&key==="current"){window.__cancelledPending=true;queueMicrotask(()=>document.querySelector('+json.dumps('button[aria-label="Close dialog"]')+')?.click())}return request}')
    click(b,'Open fragment');b.wait('window.__cancelledPending&&!document.querySelector("[role=dialog]")&&!document.querySelector(".application").inert',timeout=60)
    b.evaluate('IDBObjectStore.prototype.put=window.__put');assert slots(b)==before_cancel and b.evaluate('document.querySelectorAll(".document-entry").length')==before_documents and b.evaluate('document.querySelector("h1").textContent')=='MVP_CIRCULAR'
    # Opening is a fresh worker derivation with a retained exact source-event parent.
    click(b,'Open fragment');b.wait('!document.querySelector("[role=dialog]")&&document.querySelector("h1")?.textContent==="MVP_CIRCULAR_digest_fragment_1"',timeout=60)
    workspace=archive(b)
    with zipfile.ZipFile(io.BytesIO(workspace))as z:
        manifest=json.loads(z.read('manifest.json'));doc=next(d for d in manifest['documents']if d['id']==manifest['active_document']);event=next(e for e in manifest['history']if e['id']==doc['head_event']);assert event['parents']==[evidence['source']['binding']['event_id']]
        assert any(e['id']==evidence['source']['binding']['event_id']and e['document_id']==evidence['source']['binding']['document_id']for e in manifest['history'])
    native_product=download(b,'Save .cg2','MVP_CIRCULAR_digest_fragment_1.cg2');assert b'restriction_inspection:recipe'in native_product
    count_before_native=b.evaluate('document.querySelectorAll(".document-entry").length');upload(b,DIRECTORY/'MVP_CIRCULAR_digest_fragment_1.cg2');b.wait('document.querySelectorAll(".document-entry").length>'+str(count_before_native));click(b,'Digest');assert not b.evaluate('!!document.querySelector(".digest-results")')
    assert b.evaluate('document.querySelector(".digest-provenance")?.textContent.includes("unverified")')
    recipe=json.loads(download(b,'Download retained digest recipe','digest-retained-recipe.json'));assert recipe['fragment_index']==0
    print('PASS stage: real GenBank viewer, rotation/selection, compound annotations, independent multi-enzyme geometry, fresh product export and exact-parent native installation',flush=True)
    # Other substrate/topology cases and exact motif lengths.
    for name in ['MVP_LINEAR','MVP_EDGE','MVP_AMBIG','MVP_UNCUT','MVP_MANY','pCG_demo']:
        upload(b,paths[name]);b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name))
    assert not b.evaluate('!!document.querySelector(".mvp-example")'),'A real file name must not mark a document as an example'
    document(b,'MVP_LINEAR');assert b.evaluate('document.querySelector(".map-panel h2").textContent==="Linear map"');assert 'linear'in b.evaluate('document.querySelector(".map-stage canvas").getAttribute("aria-label")')
    linear=digest(b,['EcoRI','BamHI']);assert linear['report']['ref_top_lengths']==oracle(dna,False,['EcoRI','BamHI'])[1]
    document(b,'MVP_CIRCULAR');noti=digest(b,['NotI']);assert noti['report']['sites'][0]['recognition_length']==8
    click_expression(b,'document.querySelector(".digest-site button")');b.wait('document.querySelector(".statusbar")?.textContent.includes("101–108")')
    document(b,'MVP_EDGE');edge=digest(b,['SapI']);assert edge['report']['scan_complete']and not edge['report']['products_determined'];assert 'CutsOutsideSubstrate'in edge['report']['sites'][0]['diagnostics'];assert not b.evaluate('!!document.querySelector(".digest-fragment")')
    document(b,'MVP_AMBIG');ambig=digest(b,['EcoRI']);assert not ambig['report']['products_determined']and any(not s['certain']for s in ambig['report']['sites'])
    document(b,'MVP_UNCUT');uncut=digest(b,['EcoRI']);assert uncut['report']['products_determined']and uncut['report']['ref_top_lengths']==[47]and uncut['report']['fragments'][0]['duplex']['topology']=='Circular'and uncut['report']['fragments'][0]['duplex']['ends']is None
    document(b,'MVP_MANY');many=digest(b,['EcoRI']);assert len(many['report']['sites'])==60 and len(many['report']['fragments'])==60
    assert b.evaluate('document.querySelectorAll(".digest-site").length')==25 and b.evaluate('document.querySelectorAll(".digest-fragment").length')==25
    first=b.evaluate('document.querySelector(".digest-site").textContent');click(b,'Next digest sites');b.wait('document.querySelector(".digest-site")?.textContent!=='+json.dumps(first));click(b,'Previous digest sites')
    click(b,'Next digest fragments');b.wait('document.querySelector(".digest-fragment h4")?.textContent.includes("26")');click(b,'Previous digest fragments')
    fill(b,'Digest comparison limit','1');click(b,'Run virtual digest');b.wait('!!document.querySelector(".digest-results")');capped=json.loads(download(b,'Download digest evidence','virtual-digest.json'));assert not capped['report']['scan_complete']and not capped['report']['products_determined']and 'Comparisons'in capped['report']['limits'];fill(b,'Digest comparison limit','50000000')
    # Cancel and source switching use deterministically held outgoing jobs.
    hold_jobs(b);click(b,'Run virtual digest');b.wait('window.__heldPosts.length>0');click(b,'Cancel digest');release_jobs(b);assert not b.evaluate('!!document.querySelector(".digest-results")')
    document(b,'MVP_CIRCULAR');original=digest(b,['EcoRI','BamHI']);upload(b,paths['MVP_CIRCULAR']);b.wait('Array.from(document.querySelectorAll(".document-entry strong")).filter(e=>e.textContent==="MVP_CIRCULAR").length===2');document(b,'MVP_CIRCULAR',0)
    hold_jobs(b);click(b,'Run virtual digest');b.wait('window.__heldPosts.length>0');b.evaluate('window.__holdJobs=false');document(b,'MVP_CIRCULAR',1);release_jobs(b);assert not b.evaluate('!!document.querySelector(".digest-results")')
    duplicate=digest(b,['EcoRI','BamHI']);assert duplicate['source']['binding']['revision']==original['source']['binding']['revision']and duplicate['source']['binding']['document_id']!=original['source']['binding']['document_id'];assert duplicate['report']['batch_revision']!=original['report']['batch_revision']
    # Advanced tools remain available through visible grouped navigation.
    click(b,'PCR');b.wait('document.querySelector(".inspector-content")?.textContent.includes("Simulate PCR")');assert not b.evaluate('document.querySelector(".more-tools").open');click(b,'Features')
    b.screenshot('/tmp/collagene2-desktop-mvp-viewer.png');click(b,'Digest');click(b,'Run virtual digest');b.wait('!!document.querySelector(".digest-results")');b.evaluate('document.querySelector(".digest-results").scrollIntoView({block:"center"})');b.screenshot('/tmp/collagene2-desktop-mvp-digest.png')
    for width in [1600,1280,1200]:
        b.call('Emulation.setDeviceMetricsOverride',{'width':width,'height':1000,'deviceScaleFactor':1,'mobile':False});b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')
        metrics=b.evaluate('({viewport:innerWidth,document:document.documentElement.scrollWidth,toolbar:document.querySelector(".editor-toolbar").scrollWidth,toolbarWidth:document.querySelector(".editor-toolbar").clientWidth,inspector:document.querySelector(".inspector").scrollWidth,inspectorWidth:document.querySelector(".inspector").clientWidth})');assert metrics['document']<=width and metrics['toolbar']<=metrics['toolbarWidth']+1 and metrics['inspector']<=metrics['inspectorWidth']+1,metrics
    b.screenshot('/tmp/collagene2-desktop-mvp-laptop.png');b.call('Page.reload');b.wait('document.body?.innerText.includes("Saved in this browser")');click(b,'Digest');assert not b.evaluate('!!document.querySelector(".digest-results")')
    assert not b.exceptions and not b.console_errors,(b.exceptions,b.console_errors)
    print('PASS: focused desktop shell/visible advanced navigation; annotated GenBank import and source fidelity; circular/linear WebGL viewer; rotation and sequence/compound selections; literal multi-enzyme cuts/products; NotI8 and wrapped-site navigation; edge/ambiguity/uncut cases; bounded paging/caps/cancel; identical-content document guards; fresh GenBank products; exact-event atomic native history; reload/no restored authority;1600/1280/1200 desktop layout')


if __name__=='__main__':main()
