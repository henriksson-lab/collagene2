#!/usr/bin/env python3
"""Development-only real clipboard editing and legacy source-fidelity checks."""
import json
import struct
import time
from pathlib import Path
from browser_qa import QA_ORIGIN, Browser, mobile_enabled
from format_browser_qa import wait_file


def main():
    b=Browser();directory=Path('/tmp/collagene2-clipboard-legacy-qa');directory.mkdir(exist_ok=True)
    for path in directory.iterdir():
        if path.is_file():path.unlink()
    origin=QA_ORIGIN
    b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(directory)})
    b.call('Browser.grantPermissions',{'origin':origin,'permissions':['clipboardReadWrite','clipboardSanitizedWrite']})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':origin,'storageTypes':'indexeddb'});b.call('Page.navigate',{'url':origin});b.call('Page.bringToFront')
    b.wait("document.body?.innerText.includes('Saved in this browser')")
    b.evaluate("navigator.clipboard.writeText('SYNTHETIC_TEST_SENTINEL')")
    def upload(path):
        b.call('DOM.enable');root=b.call('DOM.getDocument')['root']['nodeId'];node=b.call('DOM.querySelector',{'nodeId':root,'selector':'input[type=file]'})['nodeId'];b.call('DOM.setFileInputFiles',{'nodeId':node,'files':[str(path)]})
    def wait_clipboard(expected):
        deadline=time.monotonic()+10
        while time.monotonic()<deadline:
            if b.evaluate('navigator.clipboard.readText().then(text=>text==='+json.dumps(expected)+')'):return
            time.sleep(.05)
        raise AssertionError('Clipboard did not contain expected synthetic sequence: '+expected)
    def copy(expected):
        b.evaluate("document.querySelector('.clipboard-tools').open=true")
        b.click('Copy selection')
        wait_clipboard(expected)
    def select(a,last):
        b.evaluate("(()=>{const canvas=document.querySelector('.sequence-canvas');canvas.focus();const r=canvas.getBoundingClientRect();for(const [i,shift] of [["+str(a)+",false],["+str(last)+",true]])canvas.dispatchEvent(new MouseEvent('click',{bubbles:true,clientX:r.left+64+i*8.5+4,clientY:r.top+30,shiftKey:shift}));})()")
    def current(expected):
        select(0,len(expected)-1);copy(expected)
    fixture=Path('fixtures/legacy/serializer-shaped-project.collagene').resolve();source=fixture.read_bytes();upload(fixture)
    b.wait("document.querySelector('h1')?.textContent==='legacy_circle'")
    b.evaluate("document.querySelectorAll('.feature-entry')[0].click()")
    b.wait("document.querySelector('.annotation-sequence')?.textContent==='GTGC'")
    copy('GTGC');b.click('Copy reverse complement');wait_clipboard('GCAC')
    assert b.evaluate("Array.from(document.querySelectorAll('button')).find(button=>button.textContent.trim()==='Cut reference selection').disabled")
    b.click('Original legacy file');assert wait_file(directory/'original-legacy.xml')==source
    assert b.evaluate("document.querySelector('.legacy-source-rows').textContent.includes('TGCAATCG')")
    # The second record's actual source rows cannot become an inferred intact duplex.
    b.evaluate("document.querySelectorAll('.document-entry')[2].click()")
    b.click('Digest');b.click('Run virtual digest')
    b.wait("!!document.querySelector('.digest-results[data-scan-complete=true][data-products-determined=false]')")
    b.wait("document.querySelector('.digest-results')?.textContent.includes('not an established intact duplex')")
    assert b.evaluate("document.querySelectorAll('.digest-fragment').length===0")
    assert not b.evaluate("Array.from(document.querySelectorAll('.digest-results button')).some(e=>e.textContent==='Open fragment'||e.textContent==='Export fragment GenBank')")
    assert b.evaluate("document.querySelector('.legacy-source-rows').textContent.includes('AC GT TAC')")
    # Explicit byte-level XDNA fixture is independent of the Rust writer/parser.
    bases=b'ACGTTAGC';comment=b'Literal XDNA source note';xdna=bytearray(112);xdna[1]=1;xdna[2]=1;xdna[28:32]=struct.pack('>I',len(bases));xdna[96:100]=struct.pack('>I',len(comment));xdna[111]=255;xdna+=bases+comment
    path=directory/'migration.xdna';path.write_bytes(xdna);upload(path)
    b.wait("document.querySelector('h1')?.textContent==='migration.xdna'")
    b.click('Original legacy file');assert wait_file(directory/'original-legacy.xdna')==xdna
    b.click('Paste sequence');b.evaluate("(()=>{const e=document.querySelector('[aria-label=\"Sequence input\"]');e.value='>clipboard_fixture\\nACGTACGT\\n';e.dispatchEvent(new Event('input',{bubbles:true}));})()")
    b.click('Open sequence');b.wait("document.querySelector('h1')?.textContent==='clipboard_fixture'")
    select(1,3);copy('CGT')
    b.click('Cut reference selection');b.wait("document.querySelector('.document-meta')?.textContent.includes('5 bp')");current('AACGT')
    b.evaluate("document.querySelector('button[title=Undo]').click()");b.wait("document.querySelector('.document-meta')?.textContent.includes('8 bp')");current('ACGTACGT')
    select(1,3);b.evaluate("navigator.clipboard.writeText('TT')");b.click('Paste replacing selection');b.wait("document.querySelector('.document-meta')?.textContent.includes('7 bp')");current('ATTACGT')
    b.evaluate("document.querySelector('button[title=Undo]').click()");b.wait("document.querySelector('.document-meta')?.textContent.includes('8 bp')")
    select(1,1);b.evaluate("navigator.clipboard.writeText('GG')");b.click('Paste before selection');b.wait("document.querySelector('.document-meta')?.textContent.includes('10 bp')");current('AGGCGTACGT')
    b.evaluate("document.querySelector('button[title=Undo]').click()");b.wait("document.querySelector('.document-meta')?.textContent.includes('8 bp')")
    select(1,3);b.evaluate("navigator.clipboard.writeText('AC?T')");b.click('Paste replacing selection');b.wait("document.querySelector('.notice')?.textContent.includes('No sequence edit was made')");current('ACGTACGT')
    # Deferred clipboard resolution is a real asynchronous race; selection changes block edits.
    select(1,3);b.evaluate("Object.defineProperty(navigator.clipboard,'readText',{configurable:true,value:()=>new Promise(resolve=>window.__resolveClipboard=resolve)})")
    b.click('Paste replacing selection');select(4,5);b.evaluate("window.__resolveClipboard('TT');delete navigator.clipboard.readText")
    b.wait("document.querySelector('.notice')?.textContent.includes('selection changed')");current('ACGTACGT')
    # A rejected clipboard write must never turn Cut into a destructive delete.
    select(1,3);b.evaluate("Object.defineProperty(navigator.clipboard,'writeText',{configurable:true,value:()=>Promise.reject(new DOMException('denied','NotAllowedError'))})")
    b.click('Cut reference selection');b.wait("document.querySelector('.notice')?.textContent.includes('Clipboard access failed')");b.evaluate("delete navigator.clipboard.writeText");current('ACGTACGT')
    # Equal molecule content does not make two independent open documents identical.
    b.click('History');b.wait("!!document.querySelector('.history-panel')")
    prior_revision=b.evaluate("document.querySelector('.history-panel > details code').textContent")
    prior_event=b.evaluate("document.querySelectorAll('.history-panel > details code')[1].textContent")
    select(1,3);b.evaluate("Object.defineProperty(navigator.clipboard,'readText',{configurable:true,value:()=>new Promise(resolve=>window.__resolveClipboard=resolve)})");b.click('Paste replacing selection')
    b.click('Paste sequence');b.evaluate("(()=>{const e=document.querySelector('[aria-label=\"Sequence input\"]');e.value='>clipboard_fixture\\nACGTACGT\\n';e.dispatchEvent(new Event('input',{bubbles:true}));})()")
    b.click('Open sequence');b.wait("document.querySelector('h1')?.textContent==='clipboard_fixture'");select(1,3)
    assert b.evaluate("document.querySelector('.history-panel > details code').textContent")==prior_revision
    assert b.evaluate("document.querySelectorAll('.history-panel > details code')[1].textContent")!=prior_event
    b.evaluate("window.__resolveClipboard('TT');delete navigator.clipboard.readText");b.wait("document.querySelector('.notice')?.textContent.includes('document or selection changed')");current('ACGTACGT')
    # Actual key events operate only when the sequence canvas owns focus.
    select(2,5);b.evaluate("navigator.clipboard.writeText('KEYBOARD_TEST_SENTINEL')")
    b.call('Input.dispatchKeyEvent',{'type':'keyDown','key':'c','code':'KeyC','windowsVirtualKeyCode':67,'modifiers':2});b.call('Input.dispatchKeyEvent',{'type':'keyUp','key':'c','code':'KeyC','windowsVirtualKeyCode':67,'modifiers':2})
    wait_clipboard('GTAC')
    b.wait("document.body?.innerText.includes('Saved in this browser')");b.call('Page.reload',{'ignoreCache':True});b.wait("document.body?.innerText.includes('Saved in this browser')");current('ACGTACGT')
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
    b.evaluate("document.querySelector('.clipboard-tools').open=true")
    if mobile_enabled():
        b.wait('document.documentElement.scrollWidth<=innerWidth')
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    b.screenshot('/tmp/collagene2-clipboard.png')
    print('PASS: real clipboard biological/reference/RC, cut/insert/replace/undo/stale race/keyboard; legacy XML/XDNA originals, aligned rows and physical guard; reload/mobile'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))
    b.socket.close()


if __name__=='__main__':main()
