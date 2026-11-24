#!/usr/bin/env python3
"""Development-only explicit FASTQ, read validation and retained oligo workflow QA."""
import json
import time
from pathlib import Path
from browser_qa import Browser


def field(b, aria, value, event='input'):
    selector=json.dumps('[aria-label="'+aria+'"]')
    b.evaluate(f"(()=>{{const e=document.querySelector({selector});e.value={json.dumps(value)};e.dispatchEvent(new Event({json.dumps(event)},{{bubbles:true}}));}})()")


def paste(b, text):
    b.click('Paste sequence')
    field(b,'Sequence input',text)
    b.click('Open sequence')


def main():
    b=Browser();directory=Path('/tmp/collagene2-reads-qa');directory.mkdir(exist_ok=True)
    for path in directory.iterdir():
        if path.is_file():path.unlink()
    b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(directory)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':'http://127.0.0.1:8080','storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':'http://127.0.0.1:8080'})
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    paste(b,'>read_reference\nTTTCGATACGGG\n')
    b.wait("document.querySelector('h1')?.textContent==='read_reference'")
    fastq='@quality_read\nCGATAC\n+\nIIIIII\n'
    paste(b,fastq)
    b.wait("!!document.querySelector('[aria-label=\"Supplied FASTQ encoding\"]')")
    assert b.evaluate("Array.from(document.querySelectorAll('button')).find(b=>b.textContent.trim()==='Import FASTQ with selected encoding').disabled")
    field(b,'Supplied FASTQ encoding','33','change')
    b.click('Import FASTQ with selected encoding')
    b.wait("document.querySelector('h1')?.textContent==='quality_read'")
    b.wait("document.body?.innerText.includes('Quality mapping matches current calls')")
    b.click('Original FASTQ record')
    b.click('Export current FASTQ')
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    end=time.monotonic()+10
    while time.monotonic()<end and not (directory/'current.fastq').exists():time.sleep(.1)
    assert (directory/'current.fastq').read_text()==fastq
    originals=[p for p in directory.iterdir() if p.name!='current.fastq']
    assert len(originals)==1 and originals[0].read_text()==fastq,originals
    b.evaluate("Array.from(document.querySelectorAll('.document-entry')).find(e=>e.textContent.includes('read_reference')).click()")
    b.click('Read validation')
    b.evaluate("(()=>{const e=document.querySelector('[aria-label=\"Validation read 2\"]');e.checked=true;e.dispatchEvent(new Event('change',{bubbles:true}));})()")
    for aria,value in [('Read quality window','3'),('Read minimum passing','3'),('Read minimum retained','1')]:field(b,aria,value)
    b.click('Validate selected reads')
    b.wait("document.querySelector('.read-consensus-preview')?.textContent==='NNNCGATACNNN'")
    b.wait("document.body?.innerText.includes('Eligible coverage 6 / 12 bases')")
    b.click('Download read validation JSON')
    b.click('Download consensus FASTA')
    end=time.monotonic()+10
    while time.monotonic()<end and not (directory/'read-validation.json').exists():time.sleep(.1)
    report=json.loads((directory/'read-validation.json').read_text())
    assert report['eligible_coverage']==6 and report['consensus']=='NNNCGATACNNN'
    assert report['reads'][0]['placement_unique_proven']
    # Settings invalidate the previously computed evidence instead of presenting stale support.
    field(b,'Read score threshold','41')
    b.wait("!document.querySelector('.read-consensus-preview')")
    b.click('Validate selected reads')
    b.wait("document.querySelector('.read-consensus-preview')?.textContent==='NNNNNNNNNNNN'")
    # Retain modified oligos as evidence; only explicitly selected plain DNA becomes a document.
    b.click('Oligo table')
    oligos='Name,Sequence,Modification,Amount,Unit\nprimer_A,CGATAC,,25,nmol\nmodified,/5Phos/CGATAC,5Phos,50,nmol\n'
    field(b,'Oligo table source',oligos)
    b.click('Read oligo columns')
    for aria,value in [('Oligo name column','0'),('Oligo ordered sequence column','1'),('Oligo modification column','2'),('Oligo amount column','3'),('Oligo unit column','4')]:field(b,aria,value,'change')
    b.click('Interpret oligo table')
    b.wait("document.body?.innerText.includes('1 convertible unmodified DNA oligos')")
    assert b.evaluate("Array.from(document.querySelectorAll('.oligo-table-wrap input[type=checkbox]')).filter(e=>e.disabled).length===1")
    assert b.evaluate("Array.from(document.querySelectorAll('button')).find(e=>e.textContent.trim()==='Next oligo rows').disabled")
    b.click('Select convertible rows on page')
    b.click('Keep collection and open selected primers')
    b.wait("document.querySelector('h1')?.textContent==='primer_A'")
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    b.call('Page.reload',{'ignoreCache':True})
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    b.wait("document.querySelector('h1')?.textContent==='primer_A'")
    b.click('History')
    b.wait("document.body?.innerText.includes('Retained oligo collections')")
    b.click('Open oligo collection')
    b.wait("!!document.querySelector('[aria-label=\"Oligo table source\"]')")
    assert b.evaluate("document.querySelector('[aria-label=\"Oligo table source\"]').value")==oligos
    b.click('Download original oligo table')
    end=time.monotonic()+10
    while time.monotonic()<end and not (directory/'original-oligos.txt').exists():time.sleep(.1)
    assert (directory/'original-oligos.txt').read_text()==oligos
    assert not b.evaluate("document.body.innerText.includes('on : click')")
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    print('PASS: explicit FASTQ encoding/evidence, worker read support and invalidation, selected oligos and collection persistence')
    b.socket.close()


if __name__=='__main__':main()
