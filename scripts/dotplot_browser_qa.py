#!/usr/bin/env python3
"""Actual WASM dot plots, independently enumerated windows and retained WebGL rendering.
Synthetic sequences only; use isolated Chrome :9224 and preview :8080.
"""
import json
import math
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from PIL import Image
from browser_qa import QA_ORIGIN, Browser, mobile_enabled
from methods_browser_qa import fill, record, check, active
from workspace_browser_qa import wait_file

DIRECTORY=Path('/tmp/collagene2-dotplot-qa')
ORIGIN=QA_ORIGIN
DNA={c:set(v) for c,v in {'A':'A','C':'C','G':'G','T':'T','U':'T','R':'AG','Y':'CT','S':'CG','W':'AT','K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}.items()}
RC=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
AA={c:{c} for c in 'ACDEFGHIKLMNPQRSTVWYOU*'}
AA.update(B=set('DN'),Z=set('EQ'),J=set('IL'),X=set('ACDEFGHIKLMNPQRSTVWYOU'))


def oracle(a,b,window,minimum,both=True,circular_a=False,circular_b=False,stride=1,protein=False):
    masks=AA if protein else DNA
    expected=[]
    for reverse in [False,True] if both else [False]:
        oriented=b.translate(RC)[::-1] if reverse else b
        for x in range(0,len(a) if circular_a else len(a)-window+1,stride):
            for y in range(0,len(b) if circular_b else len(b)-window+1,stride):
                compatible=certain=0
                for i in range(window):
                    left=masks[a[(x+i)%len(a)]];right=masks[oriented[(y+i)%len(b)]]
                    compatible+=bool(left&right);certain+=len(left)==len(right)==1 and left==right
                if compatible>=minimum:
                    expected.append((x,len(b)-1-y if reverse else y,reverse,compatible,certain,certain>=minimum))
    return sorted(expected)


def returned(r):
    return sorted((h['reference_start'],h['query_start'],h['reverse'],h['compatible_matches'],h['certain_matches'],h['certain'])for h in r['hits'])


def choose(b,name,occurrence=0):
    b.evaluate('(()=>{const e=document.querySelector("[aria-label=\\"Dot plot query\\"]");const o=Array.from(e.options).filter(o=>o.textContent.startsWith('+json.dumps(name+' [')+'))['+str(occurrence)+'];if(!o)throw Error("Missing query");e.value=o.value;e.dispatchEvent(new Event("change",{bubbles:true}));})()')


def activate(b,name,occurrence=0):
    b.evaluate('Array.from(document.querySelectorAll(".document-entry")).filter(e=>e.querySelector("strong")?.textContent==='+json.dumps(name)+')['+str(occurrence)+'].click()')
    b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name))


def paste(b,text,title):
    b.click('Paste sequence');fill(b,'Sequence input',text);b.click('Open sequence')
    b.wait('document.querySelector("h1")?.textContent==='+json.dumps(title))
    b.click('Dot plot');b.wait('!!document.querySelector(".dotplot-panel")')


def run(b):
    b.click('Calculate dot plot');b.wait('!!document.querySelector(".dotplot-result")',timeout=60)
    b.wait('document.querySelector(".dotplot-canvas")?.dataset.painted==="true"',timeout=30)


def report(b):
    path=DIRECTORY/'sequence-dotplot.json';path.unlink(missing_ok=True)
    b.click('Download dot-plot JSON');return json.loads(wait_file(path))


def inspect(b,hit):
    pages=hit['id']//25
    for _ in range(pages):b.click('Next dot-plot windows')
    b.evaluate('document.querySelector(".dotplot-hit[data-hit-id=\\"'+str(hit['id'])+'\\"]").click()')
    b.wait('document.querySelector(".dotplot-inspection")?.dataset.hitId==='+json.dumps(str(hit['id'])))
    detail=b.evaluate('(()=>{const e=document.querySelector(".dotplot-inspection");return {...e.dataset,text:e.textContent}})()')
    assert int(detail['referenceStart'])==hit['reference_start'] and int(detail['queryStart'])==hit['query_start']
    assert detail['reverse']==str(hit['reverse']).lower()


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
    fixtures=[('linear_ref','ACGTTGCA',False),('linear_query','TGCAACGT',False),('circular_ref','ACGTAC',True),('circular_query','GTACAC',True),('ambiguous_ref','ACNT',False),('clear_query','ACGT',False),('linear_ref','ACGTTGCA',False),('linear_query','TGCAACGT',False)]
    paste(b,''.join(record(*f)for f in fixtures),'linear_ref')
    assert 'linear_ref' in b.evaluate('document.querySelector("[aria-label=\\"Dot plot query\\"]").selectedOptions[0].textContent'),'Initial query is not the active source'
    choose(b,'linear_query');fill(b,'Dot plot window','3');run(b);r=report(b)
    assert r['complete'] and r['grid_complete'] and not r['sampled']
    original_revisions=(r['reference_revision'],r['query_revision'])
    assert returned(r)==oracle('ACGTTGCA','TGCAACGT',3,3),r
    assert b.evaluate('document.querySelector(".dotplot-canvas").dataset.renderer')=='webgl2'
    reverse=next(h for h in r['hits']if h['reverse']);inspect(b,reverse)
    b.click('Fit dot plot')
    b.evaluate('document.querySelector(".dotplot-canvas").scrollIntoView({block:"center"});true')
    # A real pointer click selects the nearest retained original-coordinate anchor.
    anchor=r['hits'][0]
    position=b.evaluate('(()=>{const r=document.querySelector(".dotplot-canvas").getBoundingClientRect();return{x:r.left+52+('+str(anchor['reference_start'])+'+.5)/8*(r.width-66),y:r.top+18+('+str(anchor['query_start'])+'+.5)/8*(r.height-62)}})()')
    b.call('Input.dispatchMouseEvent',{'type':'mousePressed',**position,'button':'left','clickCount':1})
    b.call('Input.dispatchMouseEvent',{'type':'mouseReleased',**position,'button':'left','clickCount':1})
    b.wait('document.querySelector(".dotplot-inspection")?.dataset.hitId==='+json.dumps(str(anchor['id'])))
    b.click('Download dot-plot SVG');svg=wait_file(DIRECTORY/'sequence-dotplot.svg').decode();root=ET.fromstring(svg)
    assert len(root.findall('.//{http://www.w3.org/2000/svg}path'))>=len(r['hits'])
    assert r['reference_revision'] in svg and r['query_revision'] in svg

    activate(b,'ambiguous_ref');choose(b,'clear_query');fill(b,'Dot plot window','4');check(b,'Dot plot both strands',False);run(b);r=report(b)
    assert returned(r)==oracle('ACNT','ACGT',4,4,both=False)
    assert len(r['hits'])==1 and not r['hits'][0]['certain'];inspect(b,r['hits'][0])
    fill(b,'Dot plot mismatches','1');assert b.evaluate('document.querySelector(".dotplot-result")===null');run(b);r=report(b)
    assert returned(r)==oracle('ACNT','ACGT',4,3,both=False) and r['hits'][0]['certain']

    activate(b,'circular_ref');choose(b,'circular_query');fill(b,'Dot plot window','4');run(b);r=report(b)
    assert returned(r)==oracle('ACGTAC','GTACAC',4,4,both=False,circular_a=True,circular_b=True)
    assert any(len(h['reference_spans'])==2 for h in r['hits']) and any(len(h['query_spans'])==2 for h in r['hits'])
    wrapped=next(h for h in r['hits']if len(h['reference_spans'])==2);inspect(b,wrapped)
    b.click('Select reference window');b.wait('document.querySelector(".sequence-panel .sr-only")?.textContent.includes("2 reference parts highlighted")')
    check(b,'Dot plot both strands',True);run(b);r=report(b)
    assert returned(r)==oracle('ACGTAC','GTACAC',4,4,both=True,circular_a=True,circular_b=True)
    reverse_wrap=next(h for h in r['hits']if h['reverse'] and len(h['query_spans'])==2)
    inspect(b,reverse_wrap)
    assert 'Query:     ACGT' in b.evaluate('document.querySelector(".dotplot-inspection pre").textContent')
    b.click('Open query window');b.wait('document.querySelector("h1")?.textContent==="circular_query"')
    b.wait('document.querySelector(".sequence-panel .sr-only")?.textContent.includes("2 reference parts highlighted")')
    activate(b,'circular_ref');b.click('Dot plot');choose(b,'circular_query');fill(b,'Dot plot window','4')
    check(b,'Dot plot both strands',False)
    fill(b,'Dot plot stride','2');run(b);r=report(b)
    assert r['sampled'] and r['grid_complete'] and not r['complete']
    assert returned(r)==oracle('ACGTAC','GTACAC',4,4,both=False,circular_a=True,circular_b=True,stride=2)
    assert 'Sampled window grid completed' in b.evaluate('document.querySelector(".dotplot-status").textContent')
    fill(b,'Dot plot stride','1');fill(b,'Dot plot work limit','1');run(b);r=report(b)
    assert not r['complete'] and not r['grid_complete'] and r['limits']
    fill(b,'Dot plot work limit','20000000');fill(b,'Dot plot hit limit','2');run(b);r=report(b)
    assert len(r['hits'])==2 and not r['complete'] and any('Hit budget' in s for s in r['limits'])
    fill(b,'Dot plot hit limit','100000');fill(b,'Dot plot index limit','1');run(b);r=report(b)
    assert not r['hits'] and not r['complete'] and any('index' in s.lower()for s in r['limits'])
    fill(b,'Dot plot index limit','1000000')

    # Identical sequence hashes do not identify a query document or reference event.
    activate(b,'linear_ref',0);choose(b,'linear_query',0);fill(b,'Dot plot window','3');check(b,'Dot plot both strands',True)
    for change in [lambda:choose(b,'linear_query',1),lambda:activate(b,'linear_ref',1)]:
        b.evaluate('window.qaPost=Worker.prototype.postMessage;Worker.prototype.postMessage=function(){};true')
        b.click('Calculate dot plot');change()
        b.evaluate('Worker.prototype.postMessage=window.qaPost;delete window.qaPost;true')
        b.wait('!Array.from(document.querySelectorAll("button")).some(e=>e.textContent==="Cancel dot plot")')
        assert b.evaluate('document.querySelector(".dotplot-result")===null')
    fill(b,'Dot plot window','3');run(b)
    current=report(b);assert (current['reference_revision'],current['query_revision'])==original_revisions
    # Deliberate WebGL unavailability uses disclosed display bins, retaining every hit.
    expected=report(b)['hits'];b.click('Features')
    b.evaluate('window.qaContext=HTMLCanvasElement.prototype.getContext;HTMLCanvasElement.prototype.getContext=function(type,...args){if(type==="webgl2"&&this.classList.contains("dotplot-canvas"))return null;return window.qaContext.call(this,type,...args)};true')
    b.click('Dot plot');choose(b,'linear_query');fill(b,'Dot plot window','3');run(b)
    assert b.evaluate('document.querySelector(".dotplot-canvas").dataset.renderer')=='canvas-aggregate'
    assert '64×64 grid' in b.evaluate('document.querySelector(".dotplot-panel").textContent')
    assert report(b)['hits']==expected
    b.evaluate('HTMLCanvasElement.prototype.getContext=window.qaContext;delete window.qaContext;true')
    b.click('Features')

    b.evaluate('window.qaTexts=[];window.qaFillText=CanvasRenderingContext2D.prototype.fillText;CanvasRenderingContext2D.prototype.fillText=function(text,...args){window.qaTexts.push(String(text));return window.qaFillText.call(this,text,...args)};true')
    paste(b,'>protein_ref [molecule=protein]\nABCD\n>protein_query [molecule=protein]\nADCD\n','protein_ref [molecule=protein]')
    assert b.evaluate('document.querySelector("[aria-label=\\"Dot plot window\\"]").value')=='4'
    assert 'protein_ref' in b.evaluate('document.querySelector("[aria-label=\\"Dot plot query\\"]").selectedOptions[0].textContent')
    assert b.evaluate('document.querySelector("[aria-label=\\"Dot plot both strands\\"]").disabled && !document.querySelector("[aria-label=\\"Dot plot both strands\\"]").checked')
    assert b.evaluate('document.querySelector(".document-entry.active small").textContent.startsWith("4 aa")')
    b.wait('window.qaTexts.some(s=>s.includes("4 aa"))')
    choose(b,'protein_query [molecule=protein]');run(b);r=report(b)
    assert returned(r)==oracle('ABCD','ADCD',4,4,both=False,protein=True) and r['alphabet']=='Protein'
    inspect(b,r['hits'][0]);b.click('Select reference window')
    b.wait('document.querySelector(".statusbar").textContent.includes("4 aa")')
    b.click('Edit selection');b.wait('document.querySelector("#dialog-title")?.textContent==="Edit selected residues"')
    assert 'Replace residues 1–4' in b.evaluate('document.querySelector("[role=dialog]").textContent')
    b.evaluate('document.querySelector("[aria-label=\\"Close dialog\\"]").click()')
    b.wait('!document.querySelector("[role=dialog]")')
    paste(b,'>rna_ref [molecule=RNA]\nAUGC\n>rna_query [molecule=RNA]\nGCAU\n','rna_ref [molecule=RNA]')
    assert b.evaluate('document.querySelector(".document-entry.active small").textContent.startsWith("4 nt")')
    assert b.evaluate('document.querySelector(".document-meta span").textContent==="4 nt"')
    b.wait('window.qaTexts.some(s=>s.includes("4 nt"))')
    b.evaluate('CanvasRenderingContext2D.prototype.fillText=window.qaFillText;delete window.qaFillText;true')
    choose(b,'rna_query [molecule=RNA]');check(b,'Dot plot both strands',True);run(b);r=report(b)
    assert len(r['hits'])==1 and r['hits'][0]['reverse'] and r['hits'][0]['query_start']==3
    inspect(b,r['hits'][0]);b.click('Select reference window')
    b.wait('document.querySelector(".statusbar").textContent.includes("4 nt")')

    # Thirty thousand bases produce an indexed, complete sparse self-plot, not a
    # quadratic matrix or scientific recomputation during camera movement.
    value=7;letters=[]
    for _ in range(30000):
        value=(1664525*value+1013904223)&0xffffffff;letters.append('ACGT'[(value>>24)&3])
    paste(b,record('large_dotplot',''.join(letters)),'large_dotplot')
    run(b);r=report(b)
    assert r['complete'] and len(r['hits'])>=30000-15+1
    assert b.evaluate('document.querySelector(".dotplot-canvas").dataset.renderer')=='webgl2'
    assert b.evaluate('document.querySelectorAll(".dotplot-hit").length')==25
    b.click('Next dot-plot windows');b.wait('document.querySelector(".dotplot-hit")?.dataset.hitId==="25"')
    b.click('Previous dot-plot windows');b.wait('document.querySelector(".dotplot-hit")?.dataset.hitId==="0"')
    b.click('Expand plot');b.evaluate('document.querySelector(".dotplot-canvas").scrollIntoView({block:"center"});true')
    b.wait('Number(document.querySelector(".dotplot-canvas").getBoundingClientRect().width)>500')
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.evaluate('window.qaUploads=0;window.qaBufferData=WebGL2RenderingContext.prototype.bufferData;WebGL2RenderingContext.prototype.bufferData=function(...args){window.qaUploads++;return window.qaBufferData.apply(this,args)};true')
    intervals=b.evaluate('new Promise(resolve=>{let previous,steps=0,intervals=[];function tick(now){if(previous!==undefined)intervals.push(now-previous);previous=now;const c=document.querySelector(".dotplot-canvas"),r=c.getBoundingClientRect();c.dispatchEvent(new WheelEvent("wheel",{deltaY:steps<15?-1:1,clientX:r.left+r.width/2,clientY:r.top+r.height/2,bubbles:true,cancelable:true}));if(++steps<30)requestAnimationFrame(tick);else requestAnimationFrame(()=>resolve(intervals));}requestAnimationFrame(tick)})')
    assert b.evaluate('window.qaUploads')==0,'Pan/zoom rebuilt retained vertex buffers'
    b.evaluate('WebGL2RenderingContext.prototype.bufferData=window.qaBufferData;delete window.qaBufferData;true')
    p95=sorted(intervals)[math.ceil(len(intervals)*.95)-1];assert p95<250,intervals
    b.click('Fit dot plot');b.screenshot('/tmp/collagene2-dotplot-desktop.png')
    rect=b.evaluate('document.querySelector(".dotplot-canvas").getBoundingClientRect().toJSON()')
    image=Image.open('/tmp/collagene2-dotplot-desktop.png').convert('RGB')
    crop=image.crop((int(rect['left']+52),int(max(0,rect['top']+18)),int(rect['right']-14),int(min(1000,rect['bottom']-44))))
    assert sum(g>r+20 and g>b+10 for r,g,b in (crop.get_flattened_data() if hasattr(crop,"get_flattened_data") else crop.getdata()))>20,'WebGL plot did not visibly paint its anchors'
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        b.wait('document.documentElement.scrollWidth<=innerWidth')
        b.screenshot('/tmp/collagene2-dotplot-mobile.png')
    b.evaluate('document.querySelector(".dotplot-canvas").dispatchEvent(new KeyboardEvent("keydown",{key:"Escape",bubbles:true}));true')
    assert not b.evaluate('document.querySelector(".dotplot-expanded")')
    assert not b.exceptions and not b.console_errors,(b.exceptions,b.console_errors)
    print('PASS: literal full-window/RC/circular/ambiguous/protein/RNA coordinates; sampling and bounds; exact source identity cancellation; exports; disclosed fallback; retained WebGL pan/zoom and visible paint; mobile. Local30kb zoom frame p95:'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"),round(p95,3),'ms')


if __name__=='__main__':main()
