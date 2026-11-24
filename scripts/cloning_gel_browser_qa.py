#!/usr/bin/env python3
"""Development Golden Gate, mutagenesis and reference-gel workflow regressions."""
import json
from pathlib import Path
from browser_qa import Browser, mobile_enabled
from format_browser_qa import wait_file


def main():
    b=Browser();directory=Path('/tmp/collagene2-cloning-gel-qa');directory.mkdir(exist_ok=True)
    for path in directory.iterdir():
        if path.is_file():path.unlink()
    b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(directory)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':'http://127.0.0.1:8080','storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':'http://127.0.0.1:8080'})
    b.wait("document.querySelector('h1')?.textContent==='pCG_demo'")
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    def field(label,value,event='input'):
        selector='[aria-label='+json.dumps(label)+']'
        b.evaluate('(()=>{const e=document.querySelector('+json.dumps(selector)+');if(!e)throw Error("Missing field");e.value='+json.dumps(value)+';e.dispatchEvent(new Event('+json.dumps(event)+',{bubbles:true}));})()')
    def checkbox(label,checked):
        selector='[aria-label='+json.dumps(label)+']'
        b.evaluate('(()=>{const e=document.querySelector('+json.dumps(selector)+');e.checked='+json.dumps(checked)+';e.dispatchEvent(new Event("change",{bubbles:true}));})()')
    def choose(name):
        b.evaluate('Array.from(document.querySelectorAll(".document-entry")).find(e=>e.querySelector("strong").textContent==='+json.dumps(name)+').click()')
        b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name))
    def option(label,name):
        selector='[aria-label='+json.dumps(label)+']'
        value=b.evaluate('Array.from(document.querySelector('+json.dumps(selector)+').options).find(e=>e.textContent==='+json.dumps(name)+').value')
        field(label,value,'change');return int(value)
    def gb(name,sequence,circular,feature=False):
        lines=[f'LOCUS       {name:<24} {len(sequence):>6} bp    DNA     {"circular" if circular else "linear  "} SYN 01-JAN-2026','FEATURES             Location/Qualifiers']
        if feature:lines+=['     CDS             1..8','                     /label="Changed source gene"','                     /translation="SOURCE_ASSERTION"','                     /codon_start=1']
        lines+=['ORIGIN','        1 '+sequence.lower(),'//'];return '\n'.join(lines)+'\n'
    source='AGTCCGATGCTACGTAGGACTTCAGACCTGTA'
    records=[('donor_A','GGTCTCAAATGCCCCGCTTAGAGACC',False,False),('donor_B','GGTCTCAGCTTTTTTAATGAGAGACC',False,False),('gate_target','AATGCCCCGCTTTTTT',True,False),('mutation_source',source,True,False),('mutation_target','AGTGGAGATGCTACGTAGGACTTCAGACCTGTA',True,False),('mutation_linear',source,False,True),('gel_A','A'*497+'GATATC'+'A'*497,False,False),('gel_B','A'*247+'GATATC'+'A'*247+'GATATC'+'A'*494,False,False)]
    path=directory/'workflow-fixtures.gb';path.write_text(''.join(gb(*record) for record in records))
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId'];node=b.call('DOM.querySelector',{'nodeId':root,'selector':'input[type=file]'})['nodeId'];b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})
    b.wait("document.querySelector('h1')?.textContent==='donor_A'")
    choose('gate_target');b.click('Golden Gate')
    checkbox('Golden Gate template 2',True);checkbox('Golden Gate template 3',True)
    b.click('Enumerate Golden Gate products');b.wait("!!document.querySelector('.gate-report')")
    b.click('Download Golden Gate report');gate=json.loads(wait_file(directory/'golden-gate-report.json'))
    assert gate['enumeration_complete'] and not gate['exhaustive_mixture_search']
    target='AATGCCCCGCTTTTTT'
    index=next(i for i,c in enumerate(gate['candidates']) if c['restriction_site_free'] and len(c['duplex']['top']['bases'])==len(target) and target in c['duplex']['top']['bases']*2)
    b.click('Review Golden Gate candidate '+str(index+1));b.click('Verify candidate')
    b.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    assert b.evaluate("document.querySelector('.verification-result').textContent.includes('one-copy')")
    b.click('Export PDF protocol');assert wait_file(directory/'cloning-protocol.pdf').startswith(b'%PDF-')
    b.click('Discard candidate');b.click('Golden Gate')
    # Component recreation requires explicitly selecting templates again.
    checkbox('Golden Gate template 2',True);checkbox('Golden Gate template 3',True)
    field('Golden Gate state limit','1');b.click('Enumerate Golden Gate products')
    b.wait("!!document.querySelector('.gate-report')")
    (directory/'golden-gate-report.json').unlink();b.click('Download Golden Gate report')
    assert not json.loads(wait_file(directory/'golden-gate-report.json'))['enumeration_complete']

    choose('mutation_linear');b.click('Mutagenesis');option('Mutagenesis template','mutation_linear')
    for label,value in [('Mutagenesis forward domain','AGACCGAT'),('Mutagenesis reverse domain','TCCAGGTC'),('Mutagenesis allowed mismatches','1'),('Mutagenesis exact 3prime bases','2')]:field(label,value)
    b.click('Simulate mismatch PCR');b.wait("document.querySelector('.mutation-pcr-result h3')?.textContent==='1 conditional PCR products'")
    b.click('Download mismatch PCR evidence');mismatch=json.loads(wait_file(directory/'mismatch-pcr.json'))
    assert mismatch['products'][0]['sequence']['bases']=='AGACCGATGCTACGTAGGACTTCAGACCTGGA'
    assert [c['source_position'] for c in mismatch['products'][0]['substitutions']]==[2,30]
    b.click('Open product reference preview');b.wait("document.querySelector('h1')?.textContent==='mutation_linear_mutant_PCR_1'")
    b.click('Features');b.wait("document.querySelector('.feature-list').textContent.includes('review required')")
    b.evaluate("Array.from(document.querySelectorAll('.feature-entry')).find(e=>e.textContent.includes('Changed source gene')).click()")
    assert b.evaluate("document.querySelector('.annotation-sequence').textContent==='AGACCGAT'")
    b.click('Export GenBank');exported=wait_file(directory/'mutation_linear_mutant_PCR_1.gb').decode()
    assert '\n                     /translation=' not in exported and 'Primer-introduced changes' in exported

    choose('mutation_target');b.click('Mutagenesis');option('Mutagenesis template','mutation_source')
    for label,value in [('Mutation edit start','4'),('Mutation edit end','6'),('Mutation inserted sequence','GGA'),('Mutation forward length','8'),('Mutation reverse length','8'),('Mutation tail split','1')]:field(label,value)
    b.click('Design inverse-PCR candidate')
    b.wait("document.querySelector('.inverse-pcr-result h3')?.textContent.includes('ConsistentUnderAssumptions')")
    b.click('Review inverse-PCR candidate');b.click('Verify candidate')
    b.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    assert b.evaluate("document.querySelector('.candidate-review').textContent.includes('Phosphorylate')")
    (directory/'cloning-protocol.pdf').unlink();b.click('Export PDF protocol')
    assert wait_file(directory/'cloning-protocol.pdf').startswith(b'%PDF-')

    b.click('Gel');a_index=option('Gel source document','gel_A')
    field('Gel digest enzymes','EcoRV');field('Gel loaded mass','120');b.click('Add gel lane')
    b_index=option('Gel source document','gel_B');b.click('Add gel lane')
    option('Gel source document','gate_target');field('Gel preparation','uncut','change');field('Gel loaded mass','');b.click('Add gel lane')
    b.click('Predict gel lanes');b.wait("!!document.querySelector('.gel-result svg')")
    b.click('Download gel SVG');svg=wait_file(directory/'predicted-gel.svg').decode();assert 'data-source-revision' in svg and '<svg' in svg
    b.click('Download gel report');gel=json.loads(wait_file(directory/'predicted-gel.json'))
    assert len(gel['lanes'])==3 and len(gel['lanes'][0]['bands'])==1
    assert [f['size_bp'] for f in gel['lanes'][0]['fragments']]==[500.,500.]
    assert [f['amount_ng'] for f in gel['lanes'][0]['fragments']]==[60.,60.]
    assert len(gel['lanes'][2]['unpositioned_fragments'])==1
    assert gel['lanes'][2]['fragments'][0]['migration'] is None
    checkbox('Diagnostic candidate '+str(a_index),True);checkbox('Diagnostic candidate '+str(b_index),True)
    field('Diagnostic digest enzymes','EcoRV');checkbox('Diagnostic double digests',False)
    b.click('Rank diagnostic digests');b.wait("!!document.querySelector('.gel-diagnostics')")
    assert b.evaluate("document.querySelector('.gel-candidate').textContent.includes('Distinguished 1 pairs')")
    b.click('Download digest ranking');ranking=json.loads(wait_file(directory/'diagnostic-digests.json'))
    assert ranking['exhaustive_search'] and ranking['candidates'][0]['distinguished_pairs']==[[0,1]]
    b.click('Use digest as gel lanes');b.wait("!document.querySelector('.gel-result')")
    field('Gel migration resolution','0.03');b.wait("!document.querySelector('.gel-diagnostics')")
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        b.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        assert b.evaluate('document.documentElement.scrollWidth<=window.innerWidth'),'Mobile outer overflow'
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    print('PASS: Golden Gate finite enumeration/review/verification/PDF; mismatch changes/reviewed annotations; inverse-PCR target-bound candidate; gel mass/comigration/circular unknown/SVG; diagnostic digest ranking; stale options/mobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__=='__main__':main()
