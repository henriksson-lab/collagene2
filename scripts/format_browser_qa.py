#!/usr/bin/env python3
"""Development-only format/source/trace/SVG checks in the isolated browser test profile."""
import base64
import json
import struct
import time
from pathlib import Path
from browser_qa import Browser, mobile_enabled


def packet(kind,data):
    return bytes([kind])+struct.pack('>I',len(data))+data


def fixtures(directory):
    source=packet(9,b'SnapGene'+struct.pack('>HHH',1,15,15))
    source+=packet(0,b'\x03AAAACCCCGGGGTTTT')
    source+=packet(10,b'<Features><Feature name="origin reverse" type="CDS" directionality="2"><Segment range="13-16"/><Segment range="1-4"/></Feature></Features>')
    source+=packet(6,b'<Notes><UseCustomMapLabel>1</UseCustomMapLabel><CustomMapLabel>snapgene_fixture</CustomMapLabel></Notes>')
    (directory/'fixture.dna').write_bytes(source)
    embl='ID   embl_fixture; SV 1; circular; genomic DNA; STD; SYN; 16 BP.\nXX\nDE   Independent browser fixture.\nFH   Key             Location/Qualifiers\nFH\nFT   CDS             complement(join(13..16,1..4))\nFT                   /label="EMBL origin"\nXX\nSQ   Sequence 16 BP; 4 A; 4 C; 4 G; 4 T; 0 other;\n     aaaaccccgg ggtttt       16\n//\n'
    (directory/'fixture.embl').write_text(embl)
    # Synthetic SCF 3 signals. Encoder uses the documented finite-difference definition.
    # Native parser tests independently use manually specified encoded values.
    channels=[[0]*14 for _ in range(4)]
    for channel,peak in enumerate([2,5,8,11]):
        channels[channel][peak-1:peak+2]=[20,100,20]
    samples=bytearray()
    for channel in channels:
        for i,value in enumerate(channel):
            samples+=struct.pack('>H',(value-2*(channel[i-1] if i else 0)+(channel[i-2] if i>1 else 0))%65536)
    calls=struct.pack('>IIII',2,5,8,11)
    calls+=b''.join(bytes(40 if i==channel else 0 for i in range(4)) for channel in range(4))
    calls+=b'ACGT'+bytes(12)
    comment=b'NAME=trace_fixture\n\0'
    header=bytearray(128);header[:4]=b'.scf'
    for offset,value in [(4,14),(8,128),(12,4),(20,4),(24,128+len(samples)),(28,len(comment)),(32,128+len(samples)+len(calls)),(40,2)]:
        struct.pack_into('>I',header,offset,value)
    header[36:40]=b'3.00'
    scf=header+samples+calls+comment
    (directory/'fixture.scf').write_bytes(scf)
    return source,embl,scf


def wait_file(path):
    end=time.monotonic()+10
    while not path.exists() and time.monotonic()<end:
        time.sleep(.1)
    assert path.exists(),path
    return path.read_bytes()


def main():
    b=Browser();directory=Path('/tmp/collagene2-format-qa');directory.mkdir(exist_ok=True)
    for path in directory.iterdir():
        if path.is_file():path.unlink()
    source,embl,scf=fixtures(directory)
    b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(directory)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':'http://127.0.0.1:8080','storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':'http://127.0.0.1:8080'})
    b.wait("document.querySelector('h1')?.textContent==='pCG_demo'")
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    b.wait("document.body.innerText.includes('Saved in this browser')")

    def open_file(path):
        b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId']
        node=b.call('DOM.querySelector',{'nodeId':root,'selector':'input[type=file]'})['nodeId']
        b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})

    open_file(directory/'fixture.dna')
    b.wait("document.querySelector('h1')?.textContent==='snapgene_fixture'")
    assert b.evaluate("document.querySelector('.document-meta').textContent.includes('Circular')")
    assert b.evaluate("document.querySelector('.feature-list').textContent.includes('origin reverse')")
    b.click('Original .dna');assert wait_file(directory/'original.dna')==source
    b.click('SVG map');svg=wait_file(directory/'sequence-map.svg').decode()
    assert '<svg' in svg and 'snapgene_fixture' in svg and 'data-revision=' in svg
    b.click('Save .cg2');saved=directory/'snapgene_fixture.cg2';wait_file(saved)
    open_file(saved);b.wait("document.querySelector('h1')?.textContent==='snapgene_fixture'")
    open_file(directory/'fixture.embl')
    b.wait("document.querySelector('h1')?.textContent==='embl_fixture'")
    b.click('EMBL');assert wait_file(directory/'sequence.embl').decode()==embl
    open_file(directory/'fixture.scf')
    b.wait("document.querySelector('h1')?.textContent==='trace_fixture'")
    b.click('Show chromatogram')
    b.wait("!!document.querySelector('.trace-canvas')")
    b.evaluate("new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))")
    assert b.evaluate("(()=>{const c=document.querySelector('.trace-canvas');const d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;return d.some((v,i)=>i%4!==3&&v<200)})()"),'Trace not painted'
    assert b.evaluate("document.querySelector('.trace-note').textContent.includes('SCF')")
    b.evaluate("document.querySelector('.trace-canvas').dispatchEvent(new KeyboardEvent('keydown',{key:'End',bubbles:true}))")
    b.wait("document.querySelector('.statusbar').textContent.includes('4–4')")
    b.screenshot('/tmp/collagene2-trace.png')
    b.click('Edit selection')
    b.evaluate("(()=>{const t=document.querySelector('textarea');t.value='A';t.dispatchEvent(new Event('input',{bubbles:true}));})()")
    b.click('Apply edit')
    b.wait("document.querySelector('.trace-error')?.textContent.includes('differs')")
    b.click('Original trace');assert wait_file(directory/'original.scf')==scf
    b.evaluate("document.querySelector('button[title=Undo]').click()")
    b.wait("!!document.querySelector('.trace-canvas')")
    b.click('Save .cg2');trace_saved=directory/'trace_fixture.cg2';wait_file(trace_saved)
    b.wait("document.body.innerText.includes('Saved in this browser')")
    b.call('Page.reload');b.wait("document.querySelector('h1')?.textContent==='trace_fixture'")
    b.click('Show chromatogram');b.wait("!!document.querySelector('.trace-canvas')")

    # Deterministic browser fetch adapter test: intercept before the external request.
    record=b'LOCUS       URL_fixture               12 bp    DNA     linear   SYN 01-JAN-2026\nDEFINITION  Synthetic source-adapter response, not an NCBI record.\nACCESSION   SYN000001\nVERSION     SYN000001.1\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaaccccgg gg\n//\n'
    b.fetch_fixture={'prefix':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi','body':record}
    b.call('Fetch.enable',{'patterns':[{'urlPattern':'https://eutils.ncbi.nlm.nih.gov/*','requestStage':'Request'}]})
    b.click('Import URL')
    b.evaluate("(()=>{const i=document.querySelector('input[aria-label=\"Sequence source URL or accession\"]');i.value='https://www.addgene.org/12345/';i.dispatchEvent(new Event('input',{bubbles:true}));})()")
    b.click('Import source')
    b.wait("document.querySelector('[role=dialog]')?.textContent.includes('requires approved API/MCP')")
    b.evaluate("(()=>{const i=document.querySelector('input[aria-label=\"Sequence source URL or accession\"]');i.value='SYN000001.1';i.dispatchEvent(new Event('input',{bubbles:true}));})()")
    b.click('Import source')
    b.wait("document.querySelector('h1')?.textContent==='URL_fixture'")
    assert b.fulfilled_urls and 'id=SYN000001.1' in b.fulfilled_urls[-1]
    assert b.evaluate("document.querySelector('.source-provenance').textContent.includes('SYN000001.1')")
    b.click('Export GenBank');gb=wait_file(directory/'URL_fixture.gb').decode()
    assert 'source:canonical_url=' in gb and 'source:content_blake3=' in gb
    b.call('Fetch.disable');b.fetch_fixture=None
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    b.call('Page.reload');b.wait("document.querySelector('h1')?.textContent==='URL_fixture'")
    if mobile_enabled():
        assert b.evaluate('document.documentElement.scrollWidth<=innerWidth'),'Mobile overflow'
    assert not b.exceptions,json.dumps(b.exceptions)
    assert not b.console_errors,json.dumps(b.console_errors)
    print('Format/source browser checks passed: SnapGene/EMBL import and evidence, SVG export, SCF chromatogram and selection, stale-call rejection/recovery, original trace download, mocked NCBI URL/provenance import and Addgene auth diagnostic.')


if __name__=='__main__':
    main()
