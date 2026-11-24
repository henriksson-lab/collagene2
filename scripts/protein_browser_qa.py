#!/usr/bin/env python3
"""Development-only protein model/UI checks against literal independent sequences."""
import json
import time
from pathlib import Path
from browser_qa import Browser, mobile_enabled


def main():
    b=Browser();directory=Path('/tmp/collagene2-protein-qa');directory.mkdir(exist_ok=True)
    b.call('Browser.setDownloadBehavior',{'behavior':'allow','downloadPath':str(directory)})
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':'http://127.0.0.1:8080','storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':'http://127.0.0.1:8080'})
    b.wait("document.querySelector('h1')?.textContent==='pCG_demo'")
    b.wait("document.body?.innerText.includes('Saved in this browser')")

    def field(label,value,event='input'):
        selector='[aria-label='+json.dumps(label)+']'
        b.evaluate("(()=>{const e=document.querySelector("+json.dumps(selector)+");e.value="+json.dumps(value)+";e.dispatchEvent(new Event("+json.dumps(event)+",{bubbles:true}));})()")

    def open_protein(name,sequence):
        b.click('Paste sequence');field('Sequence input',f'>{name} [molecule=protein]\n{sequence}')
        b.click('Open sequence');b.wait('document.querySelector("h1")?.textContent==='+json.dumps(name+' [molecule=protein]'))
        b.click('Biology');b.wait("document.querySelector('.protein-properties')")

    def calculate():
        b.click('Calculate protein properties');b.wait("document.querySelector('.protein-result')")

    def json_report():
        path=directory/'protein-properties.json';path.unlink(missing_ok=True)
        b.click('Download protein report');deadline=time.monotonic()+15
        while not path.exists() and time.monotonic()<deadline:time.sleep(.05)
        return json.loads(path.read_text())

    b.click('Biology')
    assert b.evaluate("document.querySelector('.protein-properties').textContent.includes('never guesses a protein alphabet')")
    assert not b.evaluate("Array.from(document.querySelectorAll('button')).some(e=>e.textContent==='Calculate protein properties')")
    open_protein('literal_protein','WWYYCCC*');calculate()
    assert b.evaluate("document.querySelector('.protein-average-mass').textContent==='1026.210 Da'")
    assert b.evaluate("document.querySelector('.protein-mono-mass').textContent==='1025.323401 Da'")
    assert b.evaluate("document.querySelector('.protein-epsilon').textContent.startsWith('13980 ')")
    assert b.evaluate("document.querySelector('.protein-result').dataset.complete==='true'")
    r=json_report();assert r['terminal_stop_removed'] and r['residue_count']==7 and r['source_length']==8
    assert r['options']['disulfides']=='Reduced' and r['mass']['exact_identity']
    assert abs(r['mass']['average_da']['minimum']-1026.2103)<1e-8
    assert r['charge_model']=='independent-Henderson-Hasselbalch-Epk-6.6'

    field('Protein cysteine state','maximum','change')
    assert not b.evaluate("document.querySelector('.protein-result')")
    calculate();r=json_report();assert r['disulfide_bonds']=={'minimum':1.0,'maximum':1.0}
    assert abs(r['mass']['average_da']['minimum']-1024.19442)<1e-8
    assert r['extinction_280']['epsilon_molar_inverse_cm']['minimum']==14105
    assert r['extinction_280']['absorbance_1mg_ml_1cm']
    field('Protein cysteine state','count','change');field('Intrachain disulfide bonds','2')
    b.click('Calculate protein properties');b.wait("document.querySelector('.protein-status').textContent.includes('two explicit C')")
    assert not b.evaluate("document.querySelector('.protein-result')")
    field('Intrachain disulfide bonds','1');calculate()
    field('Protein charge pH','14.1')
    assert not b.evaluate("document.querySelector('.protein-result')")
    assert b.evaluate("Array.from(document.querySelectorAll('button')).find(e=>e.textContent==='Calculate protein properties').disabled")
    field('Protein charge pH','7')
    field('Protein stop policy','reject','change')
    b.click('Calculate protein properties');b.wait("document.querySelector('.protein-status').textContent.includes('stop symbols')")
    field('Protein stop policy','terminal','change')
    field('Protein cysteine state','reduced','change')
    field('Protein modification state','unknown','change');calculate()
    r=json_report();assert r['mass'] is None and r['isoelectric_point'] is None and r['extinction_280'] is None and not r['complete']
    assert 'modification' in ' '.join(r['unavailable'])
    field('Protein modification state','unmodified','change')

    open_protein('ambiguity_fixture','BCZX');calculate()
    r=json_report();assert r['mass'] is None and r['isoelectric_point'] is None and not r['complete']
    assert set(c['symbol']for c in r['composition'])==set('BCZX')
    field('Protein ambiguity','bounds','change');calculate();r=json_report()
    assert r['mass']['average_da']['minimum']<r['mass']['average_da']['maximum']
    assert r['mass']['monoisotopic_da']['minimum']<r['mass']['monoisotopic_da']['maximum']
    assert r['extinction_280']['epsilon_molar_inverse_cm']['minimum']==0
    assert r['extinction_280']['epsilon_molar_inverse_cm']['maximum']==5500
    assert r['isoelectric_point'] is None and not r['complete']
    assert '–' in b.evaluate("document.querySelector('.protein-average-mass').textContent")

    open_protein('isobaric_fixture','J');calculate();r=json_report()
    assert r['mass']['average_da']['minimum']==r['mass']['average_da']['maximum']
    assert not r['mass']['exact_identity'] and not r['complete']
    assert abs(r['isoelectric_point']-6.1)<1e-10
    assert b.evaluate("document.querySelector('.protein-result').textContent.includes('unresolved')")

    open_protein('unsupported_residues','UO');calculate();r=json_report()
    assert r['mass'] is None and r['net_charge_at_ph'] is None and not r['complete']
    assert 'U/O' in ' '.join(r['unavailable'])
    open_protein('unbracketed_pi','R'*10000);calculate();r=json_report()
    assert r['isoelectric_point'] is None and r['net_charge_at_ph']>9999
    assert 'bracketed' in ' '.join(r['unavailable'])

    # The number limit rejects the entire operation instead of making a partial mass.
    field('Protein residue limit','3')
    b.click('Calculate protein properties');b.wait("document.querySelector('.protein-status').textContent.includes('size limit')")
    assert not b.evaluate("document.querySelector('.protein-result')")
    field('Protein residue limit','1000000')
    # Changing an option immediately after starting prevents its response from appearing.
    b.evaluate("(()=>{Array.from(document.querySelectorAll('button')).find(e=>e.textContent==='Calculate protein properties').click();const e=document.querySelector('[aria-label=\"Protein charge pH\"]');e.value='6';e.dispatchEvent(new Event('input',{bubbles:true}));})()")
    time.sleep(.4)
    assert not b.evaluate("document.querySelector('.protein-result')")
    # Unmounting a pending panel cancels its worker; it cannot write into disposed signals.
    b.evaluate("(()=>{Array.from(document.querySelectorAll('button')).find(e=>e.textContent==='Calculate protein properties').click();Array.from(document.querySelectorAll('button')).find(e=>e.textContent==='Features').click();})()")
    time.sleep(.4);assert not b.evaluate("document.querySelector('.protein-result')")
    b.click('Biology');open_protein('mobile_fixture','ACDEFGHIKLMNPQRSTVWY');calculate()
    r=json_report();assert abs(r['isoelectric_point']-7.356542538059439)<1e-10
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
        assert b.evaluate('document.documentElement.scrollWidth<=390')
    b.screenshot('/tmp/collagene2-protein-properties.png')
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    print('PASS: protein mass/pI/A280/disulfide models, stop and ambiguity assumptions, unavailable chemistry, bounds, JSON evidence, cancellation/stale results and mobile layout'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))

if __name__=='__main__':main()
