#!/usr/bin/env python3
"""Actual reference-free contig WASM, evidence retention and source lineage QA.
Synthetic reads only. Uses isolated Chrome :9224 and local preview :8080.
"""
import base64
import io
import json
import time
import zipfile
from pathlib import Path

from browser_qa import Browser, mobile_enabled
from methods_browser_qa import fill, check
from workspace_browser_qa import slots, wait_file, manifest

DIRECTORY=Path('/tmp/collagene2-contig-qa')
ORIGIN='http://127.0.0.1:8080'


def reverse(text):
    return text.translate(str.maketrans('ACGT','TGCA'))[::-1]


def genome(length, seed=17):
    value=seed;result=[]
    for _ in range(length):
        value=(1664525*value+1013904223)&0xffffffff
        result.append('ACGT'[(value>>24)&3])
    result[0]='A';result[-1]='C'
    return ''.join(result)


def paste(b,text):
    b.click('Paste sequence');fill(b,'Sequence input',text);b.click('Open sequence')


def paste_fastq(b,records):
    paste(b,''.join('@'+name+'\n'+seq+'\n+\n'+quality+'\n' for name,seq,quality in records))
    b.wait('!!document.querySelector("[aria-label=\\"Supplied FASTQ encoding\\"]")')
    fill(b,'Supplied FASTQ encoding','33','change');b.click('Import FASTQ with selected encoding')
    b.wait('!document.querySelector("[role=dialog]")')


def panel(b):
    b.click('Read assembly');b.wait('!!document.querySelector(".contig-panel")')


def select_reads(b,names):
    b.click('Clear selected reads')
    for name in names:
        b.evaluate('(()=>{const row=Array.from(document.querySelectorAll(".contig-read-choice")).find(e=>e.querySelector("strong").textContent==='+json.dumps(name)+');if(!row)throw Error("Read missing");const e=row.querySelector("input");e.checked=true;e.dispatchEvent(new Event("change",{bubbles:true}));})()')


def policy(b,unscored=False,overlap=12):
    check(b,'Contig allow unscored reads',unscored)
    fill(b,'Contig minimum overlap',str(overlap));fill(b,'Contig minimum retained bases','1')
    fill(b,'Contig minimum support','1');check(b,'Contig quality-window trimming',False)
    fill(b,'Contig maximum error percent','10');fill(b,'Contig total cell budget','128000000')
    fill(b,'Contig pair cell budget','4000000');check(b,'Contig both strands',True)


def assemble(b):
    b.click('Assemble contigs');b.wait('!!document.querySelector(".contig-result")',timeout=60)


def report(b):
    p=DIRECTORY/'contig-assembly-report.json';p.unlink(missing_ok=True)
    b.click('Download contig report');return json.loads(wait_file(p))


def close(b):
    b.evaluate('document.querySelector("button[aria-label=\\"Close dialog\\"]").click()')
    b.wait('!document.querySelector("[role=dialog]")')


def archive(b):
    p=DIRECTORY/'collagene-workspace.cg2workspace';p.unlink(missing_ok=True)
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.click('Workspace');b.wait('!!document.querySelector("[role=dialog]")');b.click('Save workspace file')
    data=wait_file(p);close(b);return data


def saved_records(data):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        result=[]
        for item in manifest(data)['attachments'].values():
            if item['kind']=='read-assembly-v1':result.append(json.loads(z.read('blobs/'+item['blob'])))
        return result


def read_record_bytes(data,kind):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        return [z.read('blobs/'+a['blob']) for a in manifest(data)['attachments'].values() if a['kind']==kind]


def current_sequence(r,target):
    matches=[c for c in r['contigs'] if len(c['placements'])==3 and c['sequence'] in (target,reverse(target))]
    assert len(matches)==1,[(c['sequence'],len(c['placements']),c['warnings'])for c in r['contigs']]
    return matches[0]


def actual_source_oracle(r,source_sequences,source_qualities,threshold=20):
    for c in r['contigs']:
        placements={p['read_index']:p for p in c['placements']}
        emitted=[]
        for column in c['columns']:
            seen=set();eligible=[]
            for o in column['observations']:
                index=o['read_index'];at=o['call_index'];seq=source_sequences[index]
                assert 0<=at<len(seq)
                base=seq[at]
                if placements[index]['reverse']:base=reverse(base)
                assert o['base']==base,(o,placements[index])
                assert o['stored_call_index']==at,'Original untrimmed fixture calls must retain exact file indices'
                assert o['raw_score']==source_qualities[index][at]
                assert o['eligible']==(source_qualities[index][at]>=threshold and base in 'ACGT')
                assert index not in seen;seen.add(index)
                if o['eligible']:eligible.append(base)
            if column['base'] is not None:emitted.append(column['base'])
            if column['state']=='SupportedBase':
                assert eligible and len(set(eligible))==1 and column['base']==eligible[0]
            elif column['state']=='SupportedGap':
                assert column['base'] is None and column['eligible_gap_reads']
            else:
                assert column['base']=='N'
        assert ''.join(emitted)==c['sequence']


def main():
    DIRECTORY.mkdir(exist_ok=True)
    for p in DIRECTORY.iterdir():
        if p.is_file():p.unlink()
    b=Browser();b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(DIRECTORY)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':ORIGIN,'storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':ORIGIN});b.wait('document.body?.innerText.includes("Saved in this browser")')
    template=genome(96)
    # Low-confidence mismatch and insertion have uniquely anchored affine gaps;
    # the inserted base differs from both neighboring template bases.
    middle=list(template[22:74]);middle[14]=next(v for v in 'ACGT' if v!=middle[14])
    insertion=next(v for v in 'ACGT' if v not in template[44:46]);middle.insert(23,insertion);middle=''.join(middle)
    qualities=[40]*len(middle);qualities[14]=0;qualities[23]=0
    sequences=[template[:50],middle,reverse(template[50:])]
    qs=[[40]*len(sequences[0]),qualities,[40]*len(sequences[2])]
    records=[(name,seq,''.join(chr(q+33)for q in quality))for name,seq,quality in zip(['quality_a','quality_b','quality_c'],sequences,qs)]
    paste_fastq(b,records);panel(b);select_reads(b,[x[0]for x in records]);policy(b)
    assemble(b);r=report(b);assert r['complete'];c=current_sequence(r,template)
    assert c['complete'] and c['layout_unique'] and c['molecule']
    actual_source_oracle(r,sequences,qs)
    assert any(p['reverse'] for p in c['placements']),'Expected reverse source placement'
    assert any(o['raw_score']==0 and not o['eligible'] for col in c['columns']for o in col['observations'])
    assert any(col['state']=='SupportedGap' for col in c['columns']),'Called insertion must remain inspectable and be removed only by eligible gap evidence'
    assert any('-'in o['reference_aligned']or'-'in o['query_aligned']for o in r['overlaps']),'Expected an actual gapped overlap'
    assert b.evaluate('document.querySelector(".contig-layout-canvas").width>0')
    # Inspect exact source calls instead of drawing a source-filled consensus.
    b.evaluate('Array.from(document.querySelectorAll("summary")).find(e=>e.textContent==="Exact layout columns").click()')
    b.evaluate('document.querySelector(".contig-column-entry").click()');b.wait('!!document.querySelector(".contig-column-detail")')
    assert 'original call' in b.evaluate('document.querySelector(".contig-column-detail").textContent')
    # Text exports reassemble and preserve the supported consensus; evidence stays
    # in JSON/native storage, not fabricated FASTA quality values.
    b.click('Download contig FASTA');fasta=wait_file(DIRECTORY/'contig-consensus.fasta').decode()
    assert ''.join(fasta.splitlines()[1:])==c['sequence']
    b.click('Download contig GenBank');gb=wait_file(DIRECTORY/'contig-consensus.gb').decode()
    assert 'contig:engine' in gb and 'unverified' in gb.lower()
    # Recording failure must leave both IndexedDB slots and current documents
    # unchanged. The successful retry records report+product atomically.
    before=manifest(archive(b));before_slots=slots(b)
    b.click('Open reconstructed contig');b.wait('!!document.querySelector(".contig-archive-panel")')
    b.evaluate("window.qaPut=IDBObjectStore.prototype.put;IDBObjectStore.prototype.put=function(v,k){if(k==='current')throw new DOMException('Synthetic assembly put failure','QuotaExceededError');return window.qaPut.call(this,v,k)};true")
    b.click('Record read assembly');b.wait('document.querySelector(".contig-archive-panel")?.textContent.includes("Assembly was not recorded")',timeout=60)
    b.evaluate('IDBObjectStore.prototype.put=window.qaPut;delete window.qaPut;true');assert slots(b)==before_slots
    b.click('Record read assembly');b.wait('!document.querySelector(".contig-archive-panel")',timeout=60)
    b.wait('document.body?.innerText.includes("Saved in this browser")')
    after_data=archive(b);after=manifest(after_data);stored=saved_records(after_data);assert len(stored)==1
    saved=stored[0];inputs={i['input_id']:i['event']for i in saved['inputs']}
    old_docs={d['id']:d for d in before['documents']};new_docs={d['id']:d for d in after['documents']}
    for doc_id,d in old_docs.items():assert new_docs[doc_id]==d,'Recording changed a source head or navigation stack'
    product=next(d for d in after['documents']if d['id']==after['active_document'])
    event=next(e for e in after['history']if e['id']==product['head_event'])
    assert set(event['parents'])=={inputs[r['reads'][p['read_index']]['id']]['event_id']for p in c['placements']}
    assert len(event['parents'])==3
    assert manifest(base64.b64decode(slots(b)['previous']))['documents']==before['documents']
    # Native portable replacement must retain the same report and exact event DAG,
    # without restoring any session verification status.
    portable=DIRECTORY/'portable-contig.cg2workspace';portable.write_bytes(after_data)
    b.click('Workspace');b.call('DOM.enable')
    root=b.call('DOM.getDocument')['root']['nodeId']
    node=b.call('DOM.querySelector',{'nodeId':root,'selector':'[role=dialog] input[type=file][accept*=cg2workspace]'})['nodeId']
    assert node
    b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(portable)]})
    b.wait('!document.querySelector("[role=dialog]")',timeout=60)
    restored=manifest(archive(b))
    assert restored['documents']==after['documents'] and restored['history']==after['history']
    b.call('Page.reload');b.wait('document.body?.innerText.includes("Saved in this browser")')
    b.click('History');b.click('View saved read assembly');b.wait('!!document.querySelector(".saved-contig-panel")')
    assert 'unverified on opening' in b.evaluate('document.querySelector(".saved-contig-authority").textContent')
    assert b.evaluate('Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Open reconstructed contig").disabled')
    b.click('Replay saved read assembly');b.wait('document.querySelector(".saved-contig-authority")?.textContent.includes("Fresh replay")',timeout=60)
    assert 'matches recorded: true' in b.evaluate('document.querySelector(".saved-contig-panel").textContent')
    close(b)
    # An advanced current source head cannot replace the historical source. A
    # source-call action forks its exact old event instead of clobbering the edit.
    b.evaluate('Array.from(document.querySelectorAll(".document-entry")).find(e=>e.querySelector("strong").textContent==="quality_a").click()')
    b.wait('document.querySelector("h1")?.textContent==="quality_a"')
    b.evaluate('document.querySelector(".sequence-canvas").focus()')
    b.call('Input.dispatchKeyEvent',{'type':'keyDown','key':'Home','code':'Home'});b.call('Input.dispatchKeyEvent',{'type':'keyUp','key':'Home','code':'Home'})
    b.click('Edit selection');b.wait('!!document.querySelector("[aria-label=\\"Replacement sequence\\"]")')
    fill(b,'Replacement sequence','T'if sequences[0][0]!='T'else'A');b.click('Apply edit')
    b.wait('!document.querySelector("[role=dialog]")');b.wait('document.body?.innerText.includes("Saved in this browser")')
    edited=manifest(archive(b));edited_source=next(d for d in edited['documents']if d['id']==inputs['read_1']['document_id'])
    assert edited_source['head_event']!=inputs['read_1']['event_id']
    b.click('History');b.click('View saved read assembly');b.wait('!!document.querySelector(".saved-contig-panel")')
    b.click('Replay saved read assembly');b.wait('document.querySelector(".saved-contig-authority")?.textContent.includes("Fresh replay")',timeout=60)
    # Expanding the read-evidence disclosure exposes its exact original call action.
    b.evaluate('Array.from(document.querySelectorAll("[role=dialog] summary")).find(e=>e.textContent==="Read trimming and evidence").click()')
    b.click('Open retained read call');b.wait('!document.querySelector(".saved-contig-panel")')
    b.wait('document.querySelector(".statusbar").textContent.includes("Selection 1–1")')
    forked=manifest(archive(b));assert len(forked['documents'])==len(edited['documents'])+1
    assert next(d for d in forked['documents']if d['id']==edited_source['id'])==edited_source
    fork=next(d for d in forked['documents']if d['id']==forked['active_document']);fork_event=next(e for e in forked['history']if e['id']==fork['head_event'])
    assert fork_event['parents']==[inputs['read_1']['event_id']]
    # Minimum-support and quality-window effects are model declarations, never
    # reference fill. Uncovered/low-support ends remain Ns.
    panel(b);select_reads(b,['quality_a','quality_b','quality_c']);policy(b);fill(b,'Contig minimum support','2')
    assemble(b);support_report=report(b)
    assert any('N'in c['sequence']for c in support_report['contigs'])
    trim_template=genome(60,51)
    trim_scores=[0]*7+[40]*44+[0]*9
    paste_fastq(b,[('trimmed_calls',trim_template,''.join(chr(v+33)for v in trim_scores))])
    panel(b);select_reads(b,['trimmed_calls']);policy(b)
    check(b,'Contig quality-window trimming',True)
    fill(b,'Contig trim window','1');fill(b,'Contig passing calls','1')
    assemble(b);tr=report(b)
    assert (tr['reads'][0]['trim_start'],tr['reads'][0]['trim_end'])==(7,51)
    tc=tr['contigs'][0];assert tc['sequence']==trim_template[7:51]
    assert tc['placements'][0]['calls'][0]['stored_call_index']==7
    actual_source_oracle(tr,[trim_template],[trim_scores])
    # Complete exact tiled plain reads, a genuinely branching repeat and explicit
    # cell exhaustion remain separately visible.
    plain=genome(96,31)
    paste(b,''.join('>'+name+'\n'+seq+'\n'for name,seq in [('plain_a',plain[:50]),('plain_b',plain[22:74]),('plain_c',plain[50:])]))
    panel(b);select_reads(b,['plain_a','plain_b','plain_c']);policy(b,True)
    assemble(b);p=current_sequence(report(b),plain);assert p['layout_unique']
    fill(b,'Contig total cell budget','1');assemble(b);limited=report(b);assert not limited['complete'] and limited['limits']
    assert 'limit' in b.evaluate('document.querySelector(".contig-status").textContent').lower()
    fill(b,'Contig total cell budget','128000000')
    # Synchronous option invalidation must cancel before a worker reply can appear.
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Assemble contigs").click();const e=document.querySelector("[aria-label=\\"Contig minimum overlap\\"]");e.value="13";e.dispatchEvent(new Event("input",{bubbles:true}));})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".contig-result")===null')
    # Equal sequence content cannot substitute for the selected source identity.
    paste(b,'>plain_equal_identity\n'+plain[:50]+'\n');panel(b)
    select_reads(b,['plain_a','plain_b','plain_c']);policy(b,True)
    b.evaluate('(()=>{Array.from(document.querySelectorAll("button")).find(e=>e.textContent==="Assemble contigs").click();for(const name of ["plain_a","plain_equal_identity"]){const row=Array.from(document.querySelectorAll(".contig-read-choice")).find(e=>e.querySelector("strong").textContent===name);const input=row.querySelector("input");input.checked=name==="plain_equal_identity";input.dispatchEvent(new Event("change",{bubbles:true}));}})()')
    time.sleep(.3);assert b.evaluate('document.querySelector(".contig-result")===null')
    repeat=[('repeat_a','ACGTACGTACGT'),('repeat_b','ACGTACGT')]
    paste(b,''.join('>'+name+'\n'+seq+'\n'for name,seq in repeat));panel(b);select_reads(b,[x[0]for x in repeat]);policy(b,True,4)
    assemble(b);rr=report(b)
    assert any(o['optimal_paths_saturated']>1 for o in rr['overlaps']) or rr['unassembled']
    assert not any(c['layout_unique']and len(c['placements'])==2 for c in rr['contigs'])
    b.click('Keep assembly evidence');b.wait('!!document.querySelector(".contig-archive-panel")');b.click('Record read assembly')
    b.wait('!document.querySelector(".contig-archive-panel")',timeout=60)
    branch_seed=genome(96,131);shared=branch_seed[30:60]
    branch=[('branch_root',branch_seed[:60]),('branch_left',shared+genome(30,831)),('branch_right',shared+genome(30,331))]
    paste(b,''.join('>'+name+'\n'+seq+'\n'for name,seq in branch))
    panel(b);select_reads(b,[x[0]for x in branch]);policy(b,True,25)
    assemble(b);br=report(b)
    assert br['warnings'] or br['unassembled'] or any(c['warnings'] or not c['layout_unique'] for c in br['contigs'])
    assert not any(c['layout_unique']and len(c['placements'])==3 for c in br['contigs']),'A branch cannot become one unique linear read layout'
    # A useful larger tiled read set uses a bounded canvas window, not one DOM cell
    # per complete alignment column and not reassembly during viewport movement.
    large=genome(1000,91)
    paste(b,''.join('>'+name+'\n'+seq+'\n'for name,seq in [('large_a',large[:600]),('large_b',large[200:800]),('large_c',large[400:])]))
    panel(b);select_reads(b,['large_a','large_b','large_c']);policy(b,True,50);assemble(b);lr=report(b);lc=current_sequence(lr,large)
    assert lc['complete'] and len(lc['columns'])==1000
    b.click('Expand contigs');b.wait('document.querySelector(".contig-panel")?.classList.contains("contig-expanded")')
    fill(b,'Contig first column','800','change')
    b.evaluate('window.qaMessages=0;window.qaPost=Worker.prototype.postMessage;Worker.prototype.postMessage=function(...args){window.qaMessages++;return window.qaPost.apply(this,args)};true')
    b.click('Next contig columns');b.click('Previous contig columns')
    time.sleep(.2);assert b.evaluate('window.qaMessages')==0
    b.evaluate('Worker.prototype.postMessage=window.qaPost;delete window.qaPost;true')
    assert b.evaluate('document.querySelectorAll(".contig-column-entry").length<=25')
    assert b.evaluate('(()=>{const row=document.querySelector(".contig-read-choice"),input=row.querySelector("input");return getComputedStyle(row).flexDirection==="row"&&input.getBoundingClientRect().width<32})()'),'Source read checkboxes must remain adjacent to labels'
    # Capture the actual source-evidence viewport, rather than the picker above it.
    b.evaluate('document.querySelector(".contig-layout-container").scrollIntoView({block:"center"})')
    b.screenshot('/tmp/collagene2-contig-desktop.png')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        assert b.evaluate('document.documentElement.scrollWidth<=innerWidth'),'Contig workflow overflows mobile viewport'
        b.evaluate('document.querySelector(".contig-layout-container").scrollIntoView({block:"center"})')
        b.screenshot('/tmp/collagene2-contig-mobile.png')
    assert not b.exceptions and not b.console_errors,(b.exceptions,b.console_errors)
    print('PASS: reference-free three-read gapped/reverse/quality consensus, exact source-call oracle, constraints/limits, historical replay and source forks, atomic report/product lineage, bounded large viewport and mobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__=='__main__':main()
