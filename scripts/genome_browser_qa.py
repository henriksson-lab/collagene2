#!/usr/bin/env python3
"""Fresh WASM large-reference QA. Literal fixtures + independent small Hamming oracle."""
import json
import time
import zipfile
import io
from pathlib import Path
from browser_qa import QA_ORIGIN, Browser
from methods_browser_qa import fill, check
from workspace_browser_qa import wait_file
from crispr_browser_qa import reverse, brute_hits, guide_oracle

DIRECTORY=Path('/tmp/collagene2-genome-qa')
ORIGIN=QA_ORIGIN
FILE='input[aria-label="Genome reference file"]'
SPACER='ACGTTAGCTAGTCAGTACGA'
FINAL='GCTATCGATACGTCAGTGCA'
C12='ACGTCAGTACGATCGTACGTAGC'


def upload(b,path,selector=FILE):
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId']
    node=b.call('DOM.querySelector',{'nodeId':root,'selector':selector})['nodeId'];assert node
    b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})


def panel(b):
    b.click('Genome references');b.wait('!!document.querySelector(".genome-panel")')
    b.wait('document.querySelector(".genome-storage-status")?.textContent.includes("separate browser cache")')


def status(b,selector,text,timeout=90):
    b.wait('document.querySelector('+json.dumps(selector)+')?.textContent.includes('+json.dumps(text)+')',timeout=timeout)


def evidence(b,label,name):
    path=DIRECTORY/name;path.unlink(missing_ok=True);b.click(label);return json.loads(wait_file(path))


def roles(b):
    return b.evaluate('Array.from(document.querySelector("[aria-label=\\"Genome source reference\\"]").options).filter(e=>e.value).map(e=>({id:e.value,text:e.textContent}))')


def source(b,role):
    fill(b,'Genome source reference',role,'change');status(b,'.genome-loading-status','Manifest identity loaded')


def ingest(b,path,name,expected=None):
    fill(b,'Genome build name',name);fill(b,'Genome expected components',json.dumps(expected or []))
    upload(b,path);status(b,'.genome-ingest-status','manifest published atomically',timeout=240)
    status(b,'.genome-loading-status','Manifest identity loaded')
    return b.evaluate('document.querySelector('+json.dumps('[aria-label="Genome source reference"]')+').value')


def interval(b,contig,start,end):
    fill(b,'Genome contig ID',contig);fill(b,'Genome interval start',str(start+1));fill(b,'Genome interval end',str(end))


def guides(b,contig,start,end,preset='spcas9'):
    interval(b,contig,start,end);fill(b,'Genome nuclease profile',preset,'change')
    b.click('Enumerate interval guides');b.wait('!!document.querySelector(".genome-guide-results")',timeout=60)
    assert b.evaluate('Array.from(document.querySelectorAll(".genome-panel button")).some(e=>e.textContent==="Next guides")')
    return evidence(b,'Export genome guides JSON','genome-guides.json')


def choose(b,spacer,rev=None):
    expression='Array.from(document.querySelectorAll(".genome-guide-entry")).find(e=>e.querySelector("strong")?.textContent==='+json.dumps(spacer)+(('&&e.querySelector("button").textContent.includes('+json.dumps('reverse' if rev else 'forward')+')')if rev is not None else '')+')'
    assert b.evaluate('!!('+expression+')'),spacer
    b.evaluate('('+expression+').querySelector("button").click()')
    b.wait('!!document.querySelector(".genome-guide-entry.selected")')


def search(b,role=None,tag=False):
    if role:fill(b,'Genome search reference',role,'change')
    check(b,'Genome tag exact intended locus',tag);fill(b,'Genome maximum mismatches','0')
    b.click('Search all imported contigs')
    b.wait('document.querySelector(".genome-search-status")?.textContent.includes("Every imported contig owner")||document.querySelector(".genome-search-status")?.textContent.includes("Every declared expected component")',timeout=300)
    return evidence(b,'Export genome search JSON','genome-search.json')


def store(b,body):
    return b.evaluate('new Promise((resolve,reject)=>{const opening=indexedDB.open("collagene-genomes",1);opening.onerror=()=>reject(opening.error);opening.onsuccess=()=>{const db=opening.result;'+body+'}})')


def counts(b):
    return store(b,'const tx=db.transaction(["catalog","staging","blocks"]);const result={};for(const name of ["catalog","staging","blocks"]){const r=tx.objectStore(name).count();r.onsuccess=()=>result[name]=r.result}tx.oncomplete=()=>{db.close();resolve(result)}')


def archive(b):
    p=DIRECTORY/'collagene-workspace.cg2workspace';p.unlink(missing_ok=True)
    b.wait('document.body?.innerText.includes("Saved in this browser")');b.click('Workspace');b.click('Save workspace file')
    data=wait_file(p);b.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()');return data


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for p in DIRECTORY.iterdir():
        if p.is_file():p.unlink()
    # 33 MiB of actual canonical DNA in one >10-million-base contig. Put both
    # strand seam footprints and a distinct final hit beyond 33 MiB.
    length=33*1024*1024+137
    dna=bytearray(b'A'*length)
    starts=[65536-11,2*65536-13,3*65536-8]
    dna[starts[0]:starts[0]+23]=(SPACER+'TGG').encode()
    dna[starts[1]:starts[1]+23]=reverse(SPACER+'TGG').encode()
    dna[starts[2]:starts[2]+27]=('TTTA'+C12).encode()
    final_start=length-73;dna[final_start:final_start+23]=(FINAL+'TGG').encode()
    circle=(SPACER+'TGG'+'AACCAACCAACCAA');circle=circle[9:]+circle[:9]
    # Small record carries definite extra matches, ambiguity, both strands and a circle.
    small=SPACER+'TGG'+'AATT'+SPACER[:5]+'N'+SPACER[6:]+'TGG'+'AATT'+reverse(SPACER+'TGG')
    large=DIRECTORY/'large.fasta'
    with large.open('wb')as f:
        f.write(b'>big\n');f.write(dna);f.write(b'\n>circle\n'+circle.encode()+b'\n')
    assert large.stat().st_size>32*1024*1024 and len(dna)>10_000_000
    many=(SPACER+'TGGAAA')*60
    small_file=DIRECTORY/'small.fasta';small_file.write_text('>small\n'+small+'\n>circle\n'+circle+'\n>many\n'+many+'\n')
    expected=[{'id':'big','source_id':'big','accession':None,'length':length,'sequence_hash':None,'topology':'Linear'},
              {'id':'circle','source_id':'circle','accession':None,'length':len(circle),'sequence_hash':None,'topology':'Circular'},
              {'id':'missing','source_id':'missing','accession':None,'length':100,'sequence_hash':None,'topology':'Linear'}]
    b=Browser();b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':ORIGIN,'storageTypes':'indexeddb'})
    # Runtime.enable may replay exceptions from the previous manual QA context.
    # Keep every exception from the new navigation, excluding only old history.
    b.exceptions.clear();b.console_errors.clear()
    b.call('Page.navigate',{'url':ORIGIN});b.wait('document.body?.innerText.includes("Saved in this browser")');panel(b)
    b.evaluate('window.__sliceSizes=[];window.__originalArrayBuffer=Blob.prototype.arrayBuffer;Blob.prototype.arrayBuffer=function(){window.__sliceSizes.push(this.size);return window.__originalArrayBuffer.call(this)}')
    large_role=ingest(b,large,'Synthetic large collection',expected)
    assert max(b.evaluate('window.__sliceSizes'))<=1024*1024
    assert counts(b)['catalog']==1 and counts(b)['staging']==0
    assert b.evaluate('Array.from(document.querySelectorAll(".genome-panel button")).some(e=>e.textContent==="Next contigs")')
    assert b.evaluate('document.querySelector(".genome-panel").scrollWidth<=document.querySelector(".genome-panel").clientWidth+1')
    print('large staged import and semantic readback passed',flush=True)
    # Source coordinates are checked independently against actual synthetic DNA,
    # including reverse calls crossing a different block seam.
    for at,preset,spacer in [(starts[0],'spcas9',SPACER),(starts[1],'spcas9',SPACER),(starts[2],'cas12a',C12)]:
        r=guides(b,'big',at-10,at+40,preset);guide_oracle({'options':r['report']['options']['guide'],'candidates':r['report']['candidates']},dna.decode())
        assert any(g['guide_sequence']==spacer and g['passes_filters']for g in r['report']['candidates'])
    r=guides(b,'circle',0,len(circle));guide_oracle({'options':r['report']['options']['guide'],'candidates':r['report']['candidates']},circle)
    assert any(g['guide_sequence']==SPACER and max(g['guide_positions'])-min(g['guide_positions'])>20 for g in r['report']['candidates'])
    r=guides(b,'big',final_start-10,final_start+30);choose(b,FINAL)
    result=search(b,large_role,True);p=result['progress']
    assert p['blocks_processed']==p['blocks_total'] and p['blocks_total']>500
    assert p['bases_verified']==length+len(circle) and p['source_integrity_verified'] and p['exhaustive_for_manifest']
    assert not p['expected_components_complete'] and not p['complete']
    assert len(result['hits'])==1 and result['hits'][0]['is_intended'] and result['hits'][0]['guide_positions'][0]==final_start
    print('fresh WASM searched every large contig owner with full hashes; missing build component stays incomplete',flush=True)
    before_interval=b.evaluate('Array.from(document.querySelectorAll('+json.dumps('[aria-label^="Genome interval"]')+')).map(e=>e.value)')
    b.click('Browse original hit locus');b.wait('document.querySelector(".genome-window")?.textContent.includes('+json.dumps(large_role)+')')
    assert b.evaluate('Array.from(document.querySelectorAll('+json.dumps('[aria-label^="Genome interval"]')+')).map(e=>e.value)')==before_interval

    # View only a bounded local canvas, including a large virtualized interval.
    interval(b,'big',0,200000);b.click('Read global window');b.wait('!!document.querySelector(".genome-sequence-canvas")')
    assert b.evaluate('document.querySelectorAll(".genome-sequence-canvas").length')==1
    b.evaluate('window.__scrollPosts=0;window.__post=Worker.prototype.postMessage;Worker.prototype.postMessage=function(...args){window.__scrollPosts++;return window.__post.apply(this,args)};const e=document.querySelector(".genome-sequence-scroll");e.scrollTop=30000;e.dispatchEvent(new Event("scroll"))')
    b.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))');assert b.evaluate('window.__scrollPosts')==0;b.evaluate('Worker.prototype.postMessage=window.__post')
    assert b.evaluate('document.querySelector(".genome-sequence-scroll").scrollHeight>30000')
    b.evaluate('document.querySelector(".genome-window").scrollIntoView({block:"center"})');b.screenshot('/tmp/collagene2-genome-desktop.png')
    # Native export verifies Blob snapshots; import same file creates a new role.
    native=DIRECTORY/'reference.cgref';native.unlink(missing_ok=True);b.click('Export reference .cgref');wait_file(native)
    assert native.read_bytes()[:8]==b'CGREF\0\x01\0'
    twin_role=ingest(b,native,'Ignored native relabel');assert twin_role!=large_role
    r=guides(b,'big',final_start-10,final_start+30);choose(b,FINAL)
    fill(b,'Genome search reference',large_role,'change');check(b,'Genome tag exact intended locus',True);b.click('Search all imported contigs')
    status(b,'.genome-search-status','exact source reference role')
    result=search(b,large_role,False);assert not any(h['is_intended']for h in result['hits'])
    print('native Blob export/import and distinct-role intended guards passed',flush=True)
    # Cancellation / identical-content role invalidation. A rendered Cancel is used.
    b.click('Search all imported contigs');b.wait('!Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Cancel genome search").disabled')
    b.click('Cancel genome search');status(b,'.genome-search-status','cancelled')
    assert not b.evaluate('document.querySelector(".genome-search-progress")?.textContent.includes("complete under all declared conditions: true")')
    b.click('Search all imported contigs');source(b,large_role);assert not b.evaluate('!!document.querySelector(".genome-hit-results")')
    # Small independent oracle: expected complete collection, extra and ambiguous hits.
    expected_small=[{'id':name,'source_id':name,'accession':None,'length':len(seq),'sequence_hash':None,'topology':topology}for name,seq,topology in [('small',small,'Linear'),('circle',circle,'Circular'),('many',many,'Linear')]]
    small_role=ingest(b,small_file,'Small oracle reference',expected_small)
    r=guides(b,'small',0,23);choose(b,SPACER,False);result=search(b,small_role,True);p=result['progress'];assert p['complete']
    profile=r['report']['options']['guide']['profile']
    actual={(h['contig_id'],tuple(h['guide_positions']),tuple(h['pam_positions']),h['reverse'],h['min_mismatches'],h['max_mismatches'],h['pam_certain'],h['definite'])for h in result['hits']}
    oracle={(name,*h)for name,seq,circ in [('small',small,False),('circle',circle,True),('many',many,False)]for h in brute_hits(seq,circ,profile,SPACER,0)}
    assert actual==oracle,(actual,oracle)
    print('complete small-reference search matches independent Hamming/IUPAC oracle',flush=True)
    assert sum(h['is_intended']for h in result['hits'])==1 and any(not h['definite']for h in result['hits'])
    b.evaluate('document.querySelector(".genome-search-progress").scrollIntoView({block:"center"})');b.screenshot('/tmp/collagene2-genome-search-desktop.png')
    # Real worker hit pages and real DOM guide pages both move past 25 entries.
    first_hit=b.evaluate('document.querySelector(".genome-hit-entry").textContent');b.click('Next hits');b.wait('document.querySelector(".genome-hit-entry")?.textContent!=='+json.dumps(first_hit));b.click('Previous hits');b.wait('document.querySelector(".genome-hit-entry")?.textContent==='+json.dumps(first_hit))
    paging=guides(b,'many',0,len(many));assert len(paging['report']['candidates'])>25
    first_guide=b.evaluate('document.querySelector(".genome-guide-entry").textContent');b.click('Next guides');b.wait('document.querySelector(".genome-guide-entry")?.textContent!=='+json.dumps(first_guide));b.click('Previous guides');b.wait('document.querySelector(".genome-guide-entry")?.textContent==='+json.dumps(first_guide))
    b.click('Browse original guide locus');b.wait('document.querySelector(".genome-window")?.textContent.includes('+json.dumps(small_role)+')')
    assert b.evaluate('document.querySelector('+json.dumps('[aria-label="Genome contig ID"]')+').value')=='many'
    # Limit interruption cannot become absence. Source role itself stays bound.
    fill(b,'Genome search comparison limit','1');b.click('Search all imported contigs');status(b,'.genome-search-status','without complete verified coverage')
    assert 'complete under all declared conditions: true'not in b.evaluate('document.querySelector(".genome-search-progress").textContent')
    capped=evidence(b,'Export genome search JSON','genome-search.json');assert not capped['progress']['exhaustive_for_manifest'] and not capped['progress']['complete']
    fill(b,'Genome search comparison limit','1000000000000');fill(b,'Genome search hit limit','1');b.click('Search all imported contigs');status(b,'.genome-search-status','without complete verified coverage')
    capped=evidence(b,'Export genome search JSON','genome-search.json');assert not capped['progress']['exhaustive_for_manifest'] and not capped['progress']['complete'] and capped['progress']['limits']
    assert 'Every imported contig owner'not in b.evaluate('document.querySelector(".genome-search-status").textContent')
    fill(b,'Genome search hit limit','4096')
    print('comparison and hit limits cannot claim exhaustive coverage',flush=True)
    # A plain FASTA collection with no expected-build declaration never upgrades
    # exhaustive imported coverage into complete build coverage.
    plain=DIRECTORY/'plain.fasta';plain.write_text('>plain\n'+SPACER+'TGG\n')
    plain_role=ingest(b,plain,'Plain imported collection')
    guides(b,'plain',0,23);choose(b,SPACER);plain_report=search(b,plain_role,True)
    assert plain_report['progress']['exhaustive_for_manifest'] and plain_report['progress']['source_integrity_verified']
    assert not plain_report['progress']['expected_components_complete'] and not plain_report['progress']['complete']
    print('plain FASTA exhaustive coverage remains a collection-only result',flush=True)
    # Forced quota/write failure must roll back its namespace and preserve old data.
    before=counts(b);b.evaluate('window.__put=IDBObjectStore.prototype.put;IDBObjectStore.prototype.put=function(...args){if(this.transaction.db.name==="collagene-genomes"&&this.name==="blocks")throw new DOMException("QA quota","QuotaExceededError");return window.__put.apply(this,args)}')
    upload(b,small_file);status(b,'.genome-ingest-status','not published');b.evaluate('IDBObjectStore.prototype.put=window.__put')
    assert counts(b)==before
    # Final publication transaction failure happens after successful readback.
    b.evaluate('window.__put=IDBObjectStore.prototype.put;IDBObjectStore.prototype.put=function(...args){if(this.transaction.db.name==="collagene-genomes"&&this.name==="catalog")throw new DOMException("QA final quota","QuotaExceededError");return window.__put.apply(this,args)}')
    upload(b,small_file);status(b,'.genome-ingest-status','not published');b.evaluate('IDBObjectStore.prototype.put=window.__put');assert counts(b)==before
    damaged=DIRECTORY/'damaged.cgref';content=bytearray(native.read_bytes());content[12]^=1;damaged.write_bytes(content)
    upload(b,damaged);status(b,'.genome-ingest-status','not published');assert counts(b)==before
    print('block-write quota, final-publication quota and corrupt native readback rollbacks passed',flush=True)
    # Small archive must carry descriptors only, and no positive search state.
    data=archive(b)
    with zipfile.ZipFile(io.BytesIO(data))as z:
        m=json.loads(z.read('manifest.json'));assert len(data)<2*1024*1024
        assert sum(entry.file_size for entry in z.infolist())<2*1024*1024
        assert 'genome-reference-link-v1'in json.dumps(m)
    # Reopening metadata never resurrects current results, even with cache intact.
    b.call('Page.reload');b.wait('document.body?.innerText.includes("Saved in this browser")');panel(b)
    assert not b.evaluate('!!document.querySelector(".genome-search-progress")')
    print('workspace descriptor archive and reload restore no current search authority',flush=True)
    source(b,large_role)
    # Deleting one actual owner block prevents a fresh full scan and .cgref export.
    store(b,'const tx=db.transaction(["catalog","blocks"],"readwrite");const r=tx.objectStore("catalog").getAll();r.onsuccess=()=>{const record=r.result.map(JSON.parse).find(r=>r.build_name==="Synthetic large collection");window.__removedCache=record.cache_id;const keys=tx.objectStore("blocks").getAllKeys();keys.onsuccess=()=>{const key=keys.result.find(k=>k.startsWith(record.cache_id+"/big/")&&k.endsWith("/0000000000000000"));tx.objectStore("blocks").delete(key)}};tx.oncomplete=()=>{db.close();resolve(true)}')
    r=guides(b,'big',final_start-10,final_start+30);choose(b,FINAL);fill(b,'Genome search reference',large_role,'change');fill(b,'Genome maximum mismatches','0');b.click('Search all imported contigs');status(b,'.genome-search-status','block is missing')
    assert 'complete under all declared conditions: true'not in b.evaluate('document.querySelector(".genome-search-progress")?.textContent||""')
    # Remove linked bulk data explicitly: metadata link survives and reattachment is exact.
    b.evaluate('(()=>{const row=Array.from(document.querySelectorAll(".genome-cache-entry")).find(e=>e.textContent.includes("Synthetic large collection"));row.querySelectorAll("button")[1].click()})()')
    status(b,'.genome-loading-status','Bulk cache missing')
    b.evaluate('Array.from(document.querySelectorAll(".genome-link")).find(e=>e.textContent.includes('+json.dumps(large_role)+')).querySelector("button").click()')
    upload(b,native);status(b,'.genome-ingest-status','manifest published atomically',timeout=240)
    assert any(r['id']==large_role for r in roles(b))
    assert not b.evaluate('!!document.querySelector(".genome-search-progress")')
    assert not b.exceptions and not b.console_errors,(b.exceptions,b.console_errors)
    print('PASS: >32 MiB actual DNA streamed into separate immutable storage; readback publication; seam/reverse/circular Cas9/Cas12a global geometry; all-owner large search; independent small oracle; explicit roles/intended identity; uncertainty/limits/cancel; native export/import; quota/hash rollback; archive/reload/missing cache/reattachment; bounded desktop canvas')


if __name__=='__main__':main()
