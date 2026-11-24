#!/usr/bin/env python3
"""Development-only actual-strand repair and oligo graph verification workflows."""
import json
import subprocess
from pathlib import Path
from browser_qa import Browser
from format_browser_qa import wait_file


def main():
    b=Browser();directory=Path('/tmp/collagene2-physical-qa');directory.mkdir(exist_ok=True)
    for path in directory.iterdir():
        if path.is_file():path.unlink()
    b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(directory)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':'http://127.0.0.1:8080','storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':'http://127.0.0.1:8080'});b.wait("document.body?.innerText.includes('Saved in this browser')")
    def field(label,value,event='input'):
        selector=json.dumps('[aria-label='+json.dumps(label)+']')
        b.evaluate(f'(()=>{{const e=document.querySelector({selector});if(!e)throw Error("Missing {label}");e.value={json.dumps(value)};e.dispatchEvent(new Event({json.dumps(event)},{{bubbles:true}}));}})()')
    def checkbox(label,checked):
        selector=json.dumps('[aria-label='+json.dumps(label)+']')
        b.evaluate(f'(()=>{{const e=document.querySelector({selector});e.checked={json.dumps(checked)};e.dispatchEvent(new Event("change",{{bubbles:true}}));}})()')
    def option(label,name):
        selector=json.dumps('[aria-label='+json.dumps(label)+']')
        value=b.evaluate(f'Array.from(document.querySelector({selector}).options).find(e=>e.textContent==={json.dumps(name)}).value')
        field(label,value,'change')
    def choose(name):
        b.evaluate('Array.from(document.querySelectorAll(".document-entry")).find(e=>e.querySelector("strong").textContent==='+json.dumps(name)+').click()')
        b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name))
    def gb(name,bases,circular,kind='DNA'):
        return f'LOCUS       {name:<24} {len(bases):>6} bp    {kind}     {"circular" if circular else "linear"}   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 {bases.lower()}\n//\n'
    records=[('fill_source','GAATTCAAAA',True,'DNA'),('fill_target','AATTCAAAAGAATT',True,'DNA'),('chew_source','CTGCAGAAAA',True,'DNA'),('chew_target','GAAAAC',True,'DNA'),('oligo_top','ACGTTAGC',False,'ss-DNA'),('oligo_bottom','GCTAACGT',False,'ss-DNA'),('anneal_target','ACGTTAGC',True,'DNA'),('wrong_target','AAAAAA',True,'DNA')]
    file=directory/'physical.gb';file.write_text(''.join(gb(*row) for row in records))
    b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId'];node=b.call('DOM.querySelector',{'nodeId':root,'selector':'input[type=file]'})['nodeId'];b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(file)]})
    b.wait("document.querySelector('h1')?.textContent==='fill_source'")
    def configure_repair(target,source,enzyme):
        choose(target);b.click('End repair / annealing');field('End repair mode','repair','change');option('End repair source',source);field('End repair enzymes',enzyme);field('End repair fragment','1');field('End repair method','blunt','change');checkbox('Circular ligation after physical operation',True)
        b.click('Compute physical intermediates');b.wait("!!document.querySelector('.end-repair-report')")
    configure_repair('fill_target','fill_source','EcoRI')
    b.click('Download end-repair evidence');report=json.loads(wait_file(directory/'end-repair.json'))
    assert report['repair']['duplex']['top']['bases']=='AATTCAAAAGAATT'
    assert report['repair']['top_added_3prime']=='AATT' and report['repair']['bottom_added_3prime']=='AATT'
    assert report['repair']['duplex']['nicks']==[]
    b.click('Review repaired-product candidate');b.click('Verify candidate')
    b.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    b.click('Export HTML protocol');assert b'End repair' in wait_file(directory/'cloning-protocol.html')
    b.click('Save candidate');assert wait_file(directory/'candidate.cg2plan')[4:6]==bytes([6,0])
    choose('wrong_target');b.click('Verify candidate');b.wait("document.querySelector('.notice')?.textContent.includes('active target differs')")
    configure_repair('chew_target','chew_source','PstI')
    b.click('Review repaired-product candidate');b.click('Verify candidate')
    b.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    def configure_anneal(phosphorylate):
        choose('anneal_target');b.click('End repair / annealing');field('End repair mode','anneal','change');option('End repair source','oligo_top');option('Annealing bottom oligo','oligo_bottom');field('Annealing minimum paired bases','4');checkbox('Top oligo ordered phosphorylated',False);checkbox('Bottom oligo ordered phosphorylated',False);checkbox('Phosphorylate physical product',phosphorylate);checkbox('Circular ligation after physical operation',True)
        b.click('Compute physical intermediates');b.wait("document.querySelector('.annealing-report')?.textContent.includes('1 annealing candidates')")
    configure_anneal(False);b.click('Download annealing evidence');anneal=json.loads(wait_file(directory/'oligo-annealing.json'))
    assert anneal['annealing']['complete'] and anneal['annealing']['unique']
    assert anneal['annealing']['candidates'][0]['duplex']['bottom']['bases']=='GCTAACGT'
    b.click('Review selected annealing candidate');b.click('Verify candidate')
    b.wait("document.querySelector('.verification-result')?.textContent.includes('Invalid')")
    assert b.evaluate("document.querySelector('.verification-result').textContent.toLowerCase().includes('phosph')")
    configure_anneal(True);b.click('Review selected annealing candidate');b.click('Verify candidate')
    b.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    # New operation evidence must be present in independently parsed embedded PDF.
    b.click('Export PDF protocol');pdf=wait_file(directory/'cloning-protocol.pdf');assert pdf.startswith(b'%PDF-')
    result=subprocess.run(['pdftotext',str(directory/'cloning-protocol.pdf'),'-'],capture_output=True,text=True,check=True).stdout
    assert 'ACGTTAGC' in result and 'GCTAACGT' in result,result
    assert 'Anneal' in result and 'Phosph' in result,result
    b.click('End repair / annealing');field('End repair mode','anneal','change');option('End repair source','oligo_top');option('Annealing bottom oligo','oligo_bottom');field('Annealing minimum paired bases','4');field('Annealing comparison budget','1');b.click('Compute physical intermediates')
    b.wait("document.querySelector('.annealing-report')?.textContent.includes('Search complete false')")
    assert b.evaluate("document.querySelector('.annealing-report').textContent.includes('unique false')")
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    b.screenshot('/tmp/collagene2-physical.png')
    print('PASS: literal fill/chew strands, independent target checks, candidatev4, annealing chemistry/phosphorylation/closure, actual PDF oligos and incomplete-search guard')
    b.socket.close()


if __name__=='__main__':main()
