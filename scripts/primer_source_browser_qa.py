#!/usr/bin/env python3
"""Development-only primer and mocked authenticated provider workflow checks."""
import json
import time
from pathlib import Path
from browser_qa import Browser
from format_browser_qa import wait_file


def main():
    b=Browser();directory=Path('/tmp/collagene2-primer-source-qa');directory.mkdir(exist_ok=True)
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

    # Hand-specified reference and reverse annotation establish actual PCR tail projection.
    record='LOCUS       primer_ui_fixture         16 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\n     misc_feature    complement(13..16)\n                     /label="Reverse source"\nORIGIN\n        1 acgaccccgg ggtaac\n//\n'
    b.click('Paste sequence');field('Sequence input',record);b.click('Open sequence')
    b.wait("document.querySelector('h1')?.textContent==='primer_ui_fixture'")
    b.click('Primers')
    field('Primer name','saved_forward');field('Primer annealing domain','ACGA');field('Primer tail','CGT')
    b.click('Analyze primer & binding')
    b.wait("document.querySelector('.primer-evidence')?.textContent.includes('1 binding hits')")
    assert b.evaluate("document.querySelector('.primer-hit').textContent.includes('0 mismatches · 0 ambiguous positions · 4 exact 3′ bases')")
    b.click('Download primer evidence')
    report=json.loads(wait_file(directory/'primer-evidence.json'))
    assert report['analysis']['ordered_sequence']=='CGTACGA'
    assert report['bindings']['hits'][0]['start']==0
    b.click('Save primer annotation')
    b.wait("Array.from(document.querySelector('[aria-label=\"Saved primer annotations\"]').options).some(o=>o.textContent==='saved_forward')")
    saved_index=b.evaluate("Array.from(document.querySelector('[aria-label=\"Saved primer annotations\"]').options).find(o=>o.textContent==='saved_forward').value")
    field('Primer annealing domain','TTTT')
    field('Saved primer annotations',saved_index,'change');b.click('Load saved primer')
    assert b.evaluate("document.querySelector('[aria-label=\"Primer annealing domain\"]').value==='ACGA'")
    assert b.evaluate("document.querySelector('[aria-label=\"Primer tail\"]').value==='CGT'")
    b.click('Analyze primer & binding');b.wait("!!document.querySelector('.primer-evidence')")
    b.click('Keep as forward primer')
    field('Primer name','saved_reverse');field('Primer annealing domain','GTTA');field('Primer tail','GAC')
    b.click('Analyze primer & binding');b.wait("document.querySelector('.primer-evidence pre')?.textContent==='GACGTTA'")
    b.click('Keep as reverse primer');b.click('Use kept pair for PCR')
    b.click('Compute products')
    b.wait("document.querySelector('h1')?.textContent==='primer_ui_fixture_PCR_1'")
    assert b.evaluate("document.querySelector('.document-meta').textContent.includes('22 bp')")
    b.click('Features')
    b.evaluate("Array.from(document.querySelectorAll('.feature-entry')).find(e=>e.querySelector('strong').textContent==='Reverse source').click()")
    b.wait("document.querySelector('.annotation-sequence')?.textContent==='GTTA'")
    assert b.evaluate("document.querySelector('.annotation-details').textContent.includes('complement(16..19)')")
    # Design candidates retain tails and can populate the actual PCR form.
    b.evaluate("Array.from(document.querySelectorAll('.document-entry')).find(e=>e.querySelector('strong').textContent==='primer_ui_fixture').click()")
    b.click('Primers')
    for label,value in [('Primer target start','1'),('Primer target end','16'),('Primer minimum length','8'),('Primer maximum length','8'),('Primer pair Tm difference','100'),('Primer minimum GC','0'),('Primer maximum GC','100'),('Primer homopolymer limit','20'),('Primer complementarity limit','20'),('Design forward tail','AA'),('Design reverse tail','TT')]:field(label,value)
    b.click('Select primer pairs')
    b.wait("document.querySelector('.primer-design h3')?.textContent==='1 ranked primer pairs'")
    b.click('Download primer pairs')
    pairs=json.loads(wait_file(directory/'primer-pairs.json'));assert pairs['pairs'][0]['product_length']==20
    field('Primer target Tm','55')
    b.wait("!document.querySelector('.primer-design')")

    # Provider fixture is explicitly mocked; neither a real token nor a live account is used.
    gb='LOCUS       ADDGENE_WEB               16 bp    DNA     circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\n     misc_feature    complement(join(13..16,1..4))\n                     /label="Provider feature"\nORIGIN\n        1 aaaaccccgg ggtttt\n//\n'
    def sequence(i):return {'sequence_id':i,'sequence_description':'Browser fixture '+str(i),'length':16,'sequence':'AAAACCCCGGGGTTTT','genbank_api_url':f'https://files.example.org/{i}.gb?signature=TRANSIENT_FIXTURE'}
    catalog={'id':12345,'name':'Mock Addgene plasmid','sequences':{'public_user_full_sequences':[sequence(41)],'public_addgene_full_sequences':[sequence(42)],'public_user_partial_sequences':[sequence(43)],'public_addgene_partial_sequences':[]}}
    b.fetch_routes=[{'prefix':'https://api.developers.addgene.org/catalog/plasmid-with-sequences/12345/','body':json.dumps(catalog).encode(),'mime':'application/json','authorization':'Token BROWSER_FIXTURE_TOKEN'}, {'prefix':'https://files.example.org/42.gb','body':gb.encode(),'authorization':None}, {'prefix':'https://files.example.org/43.gb','body':gb.replace('circular','linear  ').encode(),'authorization':None}]
    b.call('Fetch.enable',{'patterns':[{'urlPattern':'https://api.developers.addgene.org/*','requestStage':'Request'},{'urlPattern':'https://files.example.org/*','requestStage':'Request'}]})
    b.click('Import URL');field('Sequence source URL or accession','https://www.addgene.org/12345/')
    field('Addgene API token','BROWSER_FIXTURE_TOKEN');b.click('Import source')
    b.wait("document.querySelector('.addgene-choices')?.textContent.includes('Mock Addgene plasmid')")
    assert b.evaluate("document.querySelector('[aria-label=\"Addgene sequence record\"]').value===''")
    b.click('Import source')
    b.wait("document.querySelector('[role=dialog]').textContent.includes('Select an Addgene sequence record')")
    # Hold an old selection's response and change the choice while it is in flight.
    b.fetch_routes[1]['defer']=True
    field('Addgene sequence record','42','change');time.sleep(1.05);b.click('Import source')
    deadline=time.monotonic()+10
    while not getattr(b,'deferred_requests',[]) and time.monotonic()<deadline:
        b.evaluate('0');time.sleep(.05)
    assert getattr(b,'deferred_requests',[]),'Mock selection request did not start'
    field('Addgene sequence record','43','change')
    b.wait("!Array.from(document.querySelectorAll('button')).find(e=>e.textContent==='Import source').disabled")
    assert b.evaluate("document.querySelector('h1').textContent!=='ADDGENE_WEB'")
    del b.fetch_routes[1]['defer']
    field('Addgene sequence record','42','change');time.sleep(1.05);b.click('Import source')
    b.wait("document.querySelector('h1')?.textContent==='ADDGENE_WEB'")
    assert b.evaluate("document.querySelector('.source-provenance').textContent.includes('sequence 42')")
    b.click('Export GenBank');exported=wait_file(directory/'ADDGENE_WEB.gb').decode()
    assert 'source:sequence_id=42' in exported and 'source:addgene_completeness=full' in exported
    assert 'BROWSER_FIXTURE_TOKEN' not in exported and 'TRANSIENT_FIXTURE' not in exported
    b.click('Save .cg2');saved=wait_file(directory/'ADDGENE_WEB.cg2')
    assert b'BROWSER_FIXTURE_TOKEN' not in saved and b'TRANSIENT_FIXTURE' not in saved
    b.click('Import URL');field('Sequence source URL or accession','addgene:12345')
    assert b.evaluate("document.querySelector('[aria-label=\"Addgene API token\"]').value===''")
    field('Addgene API token','BROWSER_FIXTURE_TOKEN');b.click('Import source')
    b.wait("!!document.querySelector('.addgene-choices')")
    field('Addgene sequence record','43','change');time.sleep(1.05);b.click('Import source')
    b.wait("document.querySelector('.source-provenance')?.textContent.includes('partial record')")
    assert b.evaluate("document.querySelector('.document-meta').textContent.includes('Linear')")
    assert not getattr(b,'request_violations',[]),getattr(b,'request_violations',[])
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    b.call('Fetch.disable')
    print('PASS: primer metrics/binding/save/reuse/tailed-PCR annotations/design; mocked Addgene explicit records/auth routing/partial provenance/export/token lifetime')


if __name__=='__main__':main()
