#!/usr/bin/env python3
"""Development viewport measurements, reported rather than inferred performance claims.

Uses the isolated Browser test profile. Software rendering is recorded explicitly;
these local measurements do not establish hardware/browser-wide acceptance.
"""
import argparse
import json
from pathlib import Path
from browser_qa import Browser


def fixture(path, length, features):
    lines=[f'LOCUS       perf_{length:<18} {length:>9} bp    DNA     circular SYN 01-JAN-2026','FEATURES             Location/Qualifiers']
    for i in range(features):
        start=i*(length-100)//max(1,features-1)
        lines.extend([f'     misc_feature    {start+1}..{start+100}',f'                     /label="Feature_{i}"'])
    lines.append('ORIGIN')
    bases=('ACGTTGCA'*((length+7)//8))[:length]
    for i in range(0,length,60):lines.append(f'{i+1:9} '+ ' '.join(bases[j:min(j+10,i+60)] for j in range(i,min(i+60,length),10)))
    path.write_text('\n'.join(lines)+'\n//\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='/tmp/collagene2-performance.json');parser.add_argument('--frames',type=int,default=120);args=parser.parse_args()
    assert 30<=args.frames<=600
    b=Browser();directory=Path('/tmp/collagene2-performance');directory.mkdir(exist_ok=True)
    b.call('Page.addScriptToEvaluateOnNewDocument',{'source':"window.__perfLongTasks=[];try{new PerformanceObserver(l=>window.__perfLongTasks.push(...l.getEntries().map(e=>({start:e.startTime,duration:e.duration})))).observe({type:'longtask',buffered:true})}catch(e){}"})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':'http://127.0.0.1:8080','storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':'http://127.0.0.1:8080'});b.wait("document.querySelector('h1')?.textContent==='pCG_demo'");b.wait("document.body.innerText.includes('Saved in this browser')")
    result={'environment':b.evaluate("(()=>{const c=document.createElement('canvas'),g=c.getContext('webgl2'),d=g?.getExtension('WEBGL_debug_renderer_info');return{userAgent:navigator.userAgent,hardwareConcurrency:navigator.hardwareConcurrency,viewport:[innerWidth,innerHeight],renderer:d?g.getParameter(d.UNMASKED_RENDERER_WEBGL):'unavailable',devicePixelRatio}})()"),'cases':[]}
    cpu=Path('/proc/cpuinfo')
    if cpu.exists():result['environment']['cpu']=next((line.split(':',1)[1].strip() for line in cpu.read_text().splitlines() if line.startswith('model name')),'unknown')
    for length,features in [(20_000,200),(1_000_000,10_000)]:
        path=directory/f'perf_{length}.gb';fixture(path,length,features)
        b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId'];node=b.call('DOM.querySelector',{'nodeId':root,'selector':'input[type=file]'})['nodeId']
        before=b.evaluate('performance.now()');b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})
        b.wait('document.querySelector("h1")?.textContent==='+json.dumps('perf_'+str(length)))
        b.wait("document.body.innerText.includes('Saved in this browser')",timeout=60)
        b.wait("document.querySelector('.track-status')?.textContent.includes('Visible tracks')",timeout=60)
        case={'length':length,'feature_count':features,'import_persist_and_track_ready_ms':b.evaluate('performance.now()')-before,'map_rendered_features':b.evaluate("document.querySelector('.map-stage canvas').getAttribute('data-rendered-feature-count')"),'map_geometry_limited':b.evaluate("document.querySelector('.map-stage canvas').getAttribute('data-geometry-limited')"),'track_viewport_status':b.evaluate("document.querySelector('.track-status')?.textContent"),'row_height':b.evaluate("document.querySelector('.sequence-canvas').dataset.rowHeight"),'measurements':{}}
        for action in ['rotate','scroll']:
            script="""new Promise(resolve=>{
              const action=ACTION,frames=FRAMES,c=document.querySelector('.map-stage canvas'),s=document.querySelector('.sequence-scroll'),r=c.getBoundingClientRect();
              const intervals=[],dispatch=[],start=performance.now();let previous=null,i=0;
              const event=(type,x,y)=>c.dispatchEvent(new PointerEvent(type,{bubbles:true,pointerId:1,clientX:x,clientY:y}));
              if(action==='rotate')event('pointerdown',r.left+r.width*.8,r.top+r.height*.5);
              const frame=t=>{if(previous!==null)intervals.push(t-previous);previous=t;const now=performance.now();
                if(action==='rotate'){const angle=i*.05;event('pointermove',r.left+r.width*.5+Math.cos(angle)*r.width*.3,r.top+r.height*.5+Math.sin(angle)*r.height*.3)}
                else{s.scrollTop=i*Math.max(1,s.scrollHeight-s.clientHeight)/(frames-1);s.dispatchEvent(new Event('scroll',{bubbles:true}))}
                dispatch.push(performance.now()-now);if(++i<frames)requestAnimationFrame(frame);else{
                  if(action==='rotate')event('pointerup',r.left+r.width*.8,r.top+r.height*.5);
                  const stats=a=>{a.sort((x,y)=>x-y);return{p50:a[Math.floor(a.length*.5)],p95:a[Math.floor(a.length*.95)],max:a[a.length-1]}};
                  requestAnimationFrame(()=>resolve({frames,frame_interval_ms:stats(intervals),event_dispatch_ms:stats(dispatch),long_tasks:window.__perfLongTasks.filter(e=>e.start>=start),canvas:[document.querySelector('.sequence-canvas').width,document.querySelector('.sequence-canvas').height],heap_js:performance.memory?.usedJSHeapSize??null}));
                }};requestAnimationFrame(frame);
            })""".replace('ACTION',json.dumps(action)).replace('FRAMES',str(args.frames))
            case['measurements'][action]=b.evaluate(script)
        if features>500:
            b.click('Next annotations')
            b.wait("document.querySelector('.feature-list strong')?.textContent==='Feature_500'")
            assert b.evaluate("document.querySelectorAll('.feature-entry').length===500")
            b.click('Previous annotations')
            b.wait("document.querySelector('.feature-list strong')?.textContent==='Feature_0'")
            case['annotation_pagination_checked']=True
        result['cases'].append(case)
    assert not b.exceptions,b.exceptions
    Path(args.output).write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
