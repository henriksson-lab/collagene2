"""Real WASM viewport/cache and unified-search checks with independent literal DNA.
Development-only Chrome CDP helper; no external scientific runtime dependencies.
"""
import json
import time
from pathlib import Path
from browser_qa import QA_ORIGIN, Browser, mobile_enabled


def set_value(browser, selector, value, event="input"):
    browser.evaluate("(()=>{const e=document.querySelector("+json.dumps(selector)+");e.value="+json.dumps(value)+";e.dispatchEvent(new Event("+json.dumps(event)+",{bubbles:true}));})()")


def toggle(browser, label, checked=True):
    browser.evaluate("(()=>{const e=document.querySelector('input[aria-label="+json.dumps(label)+"]');if(e.checked!=="+json.dumps(checked)+")e.click();})()")


def paste(browser, text, name):
    browser.click("Paste sequence")
    set_value(browser, 'textarea[aria-label="Sequence input"]', text)
    browser.click("Open sequence")
    browser.wait("document.querySelector('h1')?.textContent==="+json.dumps(name))
    browser.wait("document.querySelector('.track-status')?.textContent.includes('Visible tracks')")


def find(browser, query, mode):
    if not browser.evaluate("!!document.querySelector('select[aria-label=\"Search source\"]')"):
        browser.evaluate("document.querySelector('button[title=\"Search options and results\"]').click()")
    set_value(browser, 'select[aria-label="Search source"]', mode, 'change')
    set_value(browser, '.search input', query)
    browser.click('Find')
    browser.wait("!!document.querySelector('.search-summary')")
    return browser.evaluate("Array.from(document.querySelectorAll('.search-hit')).map(e=>e.textContent)")


def report(browser, downloads):
    browser.wait("document.querySelector('.track-status')?.textContent.includes('Visible tracks')")
    browser.evaluate("document.querySelector('.viewport-tracks').open=true")
    path=downloads/'viewport-tracks.json'
    path.unlink(missing_ok=True)
    browser.click('Download viewport report')
    deadline=time.monotonic()+15
    while not path.exists() and time.monotonic()<deadline:
        time.sleep(.1)
    return json.loads(path.read_text())


def main():
    b=Browser()
    d=Path('/tmp/collagene2-search-tracks-downloads');d.mkdir(exist_ok=True)
    b.call('Browser.setDownloadBehavior', {'behavior':'allow','downloadPath':str(d)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':QA_ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.addScriptToEvaluateOnNewDocument',{'source':'''window.__wire=[];(()=>{let next=0;const ids=new WeakMap(),original=Worker.prototype.postMessage;Worker.prototype.postMessage=function(data,...args){if(!ids.has(this))ids.set(this,++next);if(data instanceof Uint8Array)window.__wire.push({worker:ids.get(this),bytes:data.length,version:data[4]+256*data[5]});return original.call(this,data,...args);};})();'''})
    b.call('Page.navigate',{'url':QA_ORIGIN})
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    sequence='AAAAGAATTCTTTTATGAAACCCGGGTAA'
    gb=f'''LOCUS       track_fixture {len(sequence)} bp DNA linear
DEFINITION  Independent browser track fixture.
FEATURES             Location/Qualifiers
     gene            1..9
                     /label="alpha_tag"
                     /note="needle qualifier"
     primer_bind     complement(2..7)
                     /label="primer_tag"
     CDS             15..29
                     /label="known_cds"
                     /codon_start=1
                     /transl_table=11
     misc_feature    join(1..3,10..12)
                     /label="compound_tag"
ORIGIN
        1 {sequence.lower()}
//
'''
    paste(b,gb,'track_fixture')
    assert len(find(b,'needle','Annotations'))==1
    assert b.evaluate("document.querySelector('.statusbar').textContent.includes('1–9')")
    assert len(find(b,'primer_tag','Primers'))==1
    assert 'reverse' in b.evaluate("document.querySelector('.search-hit').textContent")
    assert len(find(b,'EcoRI','Enzymes'))==1
    assert 'cuts: top after base 5, bottom after base 9' in b.evaluate("document.querySelector('.search-hit').textContent")
    hits=find(b,'GAATTC','Sequence');assert len(hits)==1
    assert b.evaluate("document.querySelector('.statusbar').textContent.includes('5–10')")
    b.evaluate("document.querySelector('button[aria-label=\"Close search results\"]').click()")
    b.click('Sequence tracks')
    for label in ['Enzyme site tracks','Enzyme cut tracks','Translation tracks','CDS tracks']:
        toggle(b,label)
    set_value(b,'input[aria-label="Track minimum ORF length"]','1','change');toggle(b,'ORF tracks')
    r=report(b,d)
    assert r['complete'],r['warnings']
    assert any(i['kind']=='Primer' and i['source_feature']==1 and i['reverse'] for i in r['items'])
    eco=next(i for i in r['items'] if i['enzyme']=='EcoRI')
    assert eco['parts']==[{'start':4,'end':10}],eco
    assert eco['cuts']=={'top':5,'bottom':9},eco
    cds=next(i for i in r['items'] if i['kind']=='Cds')
    assert cds['protein']=='MKPG',cds
    assert cds['residues'][0]['codon_parts']==[{'start':14,'end':17}],cds
    assert any(i['kind']=='Orf' for i in r['items'])
    assert b.evaluate("Number(document.querySelector('.sequence-canvas').dataset.paintedTracks)>0")
    assert b.evaluate("(()=>{const c=document.querySelector('.sequence-canvas');return c.getContext('2d').getImageData(0,0,c.width,c.height).data.some(v=>v!==255&&v!==0)})()")
    # A real canvas click on forward-frame residue 1 maps AAA back to bases1..3.
    b.evaluate("(()=>{const c=document.querySelector('.sequence-canvas'),r=c.getBoundingClientRect();c.dispatchEvent(new MouseEvent('click',{bubbles:true,clientX:r.left+68,clientY:r.top+202}));})()")
    b.wait("document.querySelector('.track-selection')?.textContent.includes('residue 1 (K)')")
    assert b.evaluate("document.querySelector('.statusbar').textContent.includes('1–3')")
    b.screenshot('/tmp/collagene2-sequence-tracks.png')
    # Changing settings clears old report; the same worker serves tiny follow-up requests.
    set_value(b,'select[aria-label="Track direction"]','Reverse','change')
    r=report(b,d);assert r['options']['orientation']=='Reverse'
    assert all(i['reverse'] for i in r['items'] if i['kind']=='Translation')
    # Restore basic tracks before scrolling a larger reference. Source install transfers once.
    for label in ['Enzyme site tracks','Enzyme cut tracks','Translation tracks','CDS tracks','ORF tracks']:
        toggle(b,label,False)
    paste(b,'>large_reference\n'+'ACGT'*4000,'large_reference')
    b.wait("document.querySelector('.track-status')?.textContent==='Visible tracks complete'")
    initial=b.evaluate('window.__wire.slice()')
    groups={}
    for packet in initial:groups.setdefault(packet['worker'],[]).append(packet)
    candidates=[packets for packets in groups.values() if len(packets)>1 and packets[0]['bytes']>16000 and packets[1]['bytes']<512]
    assert candidates,initial
    install=candidates[-1][0];worker=install['worker']
    assert install['bytes']>16000,install
    for top in [2500,5000,8000,12000]:
        b.evaluate(f"(()=>{{const e=document.querySelector('.sequence-scroll');e.scrollTop={top};e.dispatchEvent(new Event('scroll'));}})()")
        time.sleep(.35)
    b.wait("document.querySelector('.track-status')?.textContent==='Visible tracks complete'")
    follow=b.evaluate('window.__wire.slice('+str(len(initial))+')')
    cached=[x for x in follow if x['worker']==worker]
    assert cached and all(x['bytes']<512 for x in cached),cached
    assert all(x['version']==22 for x in initial+follow),initial+follow
    # Rapid scroll/option changes must settle on the latest viewport and frame settings.
    b.evaluate("(()=>{const e=document.querySelector('.sequence-scroll');for(let i=0;i<30;i++){e.scrollTop=i*1000;e.dispatchEvent(new Event('scroll'));}})()")
    toggle(b,'Translation tracks')
    set_value(b,'select[aria-label="Track direction"]','Forward','change')
    r=report(b,d);assert r['options']['orientation']=='Forward' and r['options']['translations']
    assert all(i['reverse'] is False for i in r['items'] if i['kind']=='Translation')
    assert r['window']['start']>0
    # Bounded all-mode search exposes pages instead of silently dropping display rows.
    hits=find(b,'ACG','Sequence');assert len(hits)==50
    assert 'incomplete' in b.evaluate("document.querySelector('.search-summary').textContent")
    first=b.evaluate("document.querySelector('.search-hit').dataset.hitId")
    b.click('Next results')
    assert b.evaluate("document.querySelector('.search-hit').dataset.hitId")!=first
    b.evaluate("document.querySelector('button[aria-label=\"Close search results\"]').click()")
    # Circular origin hit keeps two highlights and cannot become a fabricated contiguous edit.
    paste(b,'LOCUS       circle 4 bp DNA circular\nFEATURES             Location/Qualifiers\nORIGIN\n        1 tcaa\n//\n','circle')
    hits=find(b,'ATT','Sequence');assert len(hits)==1 and 'all parts highlighted' in hits[0],hits
    assert b.evaluate("document.querySelector('.statusbar').textContent.includes('Click a base')")
    b.evaluate("document.querySelector('button[aria-label=\"Close search results\"]').click()")
    # More than50 overlapping features remain reachable in the accessible track list.
    many='LOCUS       many_features 12 bp DNA linear\nFEATURES             Location/Qualifiers\n'
    many+=''.join(f'     gene            1..3\n                     /label="overlap_{i}"\n' for i in range(61))
    many+='ORIGIN\n        1 aaacccgggttt\n//\n'
    paste(b,many,'many_features')
    b.evaluate("document.querySelector('.viewport-tracks').open=true")
    b.wait("document.querySelectorAll('.viewport-item').length===50")
    first_item=b.evaluate("document.querySelector('.viewport-item').textContent")
    b.click('Next track items')
    b.wait("document.querySelectorAll('.viewport-item').length<50")
    assert b.evaluate("document.querySelector('.viewport-item').textContent")!=first_item
    assert not b.evaluate("document.querySelector('.viewport-tracks').textContent.includes('on : click')")
    # Protein search must not pretend to search a complementary strand.
    paste(b,'>protein [molecule=protein]\nMANA','protein [molecule=protein]')
    assert len(find(b,'ANA','Sequence'))==1
    b.evaluate("document.querySelector('button[aria-label=\"Close search results\"]').click()")
    # Mobile controls/popover stay within the viewport.
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
    b.evaluate("document.querySelector('button[title=\"Search options and results\"]').click()")
    b.wait("document.querySelector('.search-results')")
    if mobile_enabled():
        assert b.evaluate('document.documentElement.scrollWidth<=390')
        assert b.evaluate("document.querySelector('.search-results').getBoundingClientRect().right<=390")
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    print('PASS: unified search/compound and protein navigation, actual cuts/CDS/codon/ORF tracks, cached binary viewport requests, stale option rejection, pagination and mobile controls'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))

if __name__=='__main__':main()
