"""Independent reciprocal BP products through the real WASM UI.

Uses the dedicated QA browser. COLLAGENE_RECOMBINATION_QA_URL can point to a
separate staging preview while the user reviews the published build on port8080.
"""
import json
import os
import subprocess
from pathlib import Path

from browser_qa import QA_ORIGIN, Browser
from methods_browser_qa import active, check, fill
from workspace_browser_qa import wait_file

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = Path("/tmp/collagene2-recombination-qa")
BASE_URL = os.environ.get("COLLAGENE_RECOMBINATION_QA_URL", QA_ORIGIN)


def choose_source(browser, label, name):
    selector = "select[aria-label=" + json.dumps(label) + "]"
    browser.evaluate("""(()=>{
        const e=document.querySelector(%s);
        const option=Array.from(e.options).find(o=>o.textContent.startsWith(%s+' · '));
        if(!option)throw Error('Missing source option');
        e.value=option.value;e.dispatchEvent(new Event('change',{bubbles:true}));
    })()""" % (json.dumps(selector), json.dumps(name)))


def panel(browser):
    browser.click("Recombination")
    browser.wait("!!document.querySelector('.recombination-panel')")


def configure(browser):
    choose_source(browser, "Recombination first source", "BP_insert")
    choose_source(browser, "Recombination second source", "BP_donor")
    check(browser, "Use manufacturer recombination context profile")


def compute(browser):
    browser.click("Reconstruct reciprocal components")
    browser.wait("document.querySelectorAll('input[name=recombination-component]').length===2")


def evidence(browser):
    path = DIRECTORY / "recombination-reference-evidence.json"
    path.unlink(missing_ok=True)
    browser.click("Download reference evidence JSON")
    return json.loads(wait_file(path))


def select_length(browser, length):
    browser.evaluate("""(()=>{
        const e=Array.from(document.querySelectorAll('input[name=recombination-component]'))
          .find(e=>e.getAttribute('aria-label').includes(%s));
        if(!e)throw Error('Missing component');e.click();
    })()""" % json.dumps(" · %s bases · " % length))


def verify_selected(browser, expected):
    browser.click("Review selected component in Plan")
    browser.wait("!!document.querySelector('.candidate-review')")
    browser.click("Verify candidate")
    browser.wait("document.querySelector('.verification-result strong')?.textContent===" + json.dumps(expected))


def main():
    DIRECTORY.mkdir(exist_ok=True)
    browser = Browser()
    browser.call("Browser.setDownloadBehavior", {"behavior": "allow", "downloadPath": str(DIRECTORY)})
    browser.call("Emulation.setDeviceMetricsOverride", {"width": 1600, "height": 1000, "deviceScaleFactor": 1, "mobile": False})
    browser.call("Storage.clearDataForOrigin", {"origin": BASE_URL, "storageTypes": "indexeddb"})
    browser.call("Page.navigate", {"url": BASE_URL})
    browser.wait("document.body?.innerText.includes('Saved in this browser')")
    browser.click("Paste sequence")
    fill(browser, "Sequence input", (ROOT / "fixtures/recombination-bp-demo.gb").read_text())
    browser.click("Open sequence")
    browser.wait("document.querySelector('h1')?.textContent==='BP_insert'")
    active(browser, "BP_target")
    panel(browser)
    assert browser.evaluate("Array.from(document.querySelectorAll('button')).find(b=>b.textContent.trim()==='Reconstruct reciprocal components').disabled")
    configure(browser)
    assert browser.evaluate("getComputedStyle(document.querySelector('.recombination-panel .checkbox-label')).flexDirection==='row'")
    compute(browser)
    report = evidence(browser)
    reference = report["reference"]
    assert reference["complete"] and len(reference["candidates"]) == 1
    products = reference["candidates"][0]["products"]
    assert sorted(len(p["molecule"]["sequence"]["bases"]) for p in products) == [2267, 2553]
    assert sum(len(p["molecule"]["sequence"]["bases"]) for p in products) == 59 + 4761
    # This independent oracle uses the original manufacturer snapshot and the
    # literal map switch boundaries, without reading reported junctions/maps.
    donor = "".join((ROOT / "crates/recombination/tests/fixtures/pdonr221_seq.txt").read_text().split()).upper()
    insert = "ACAAGTTTGTACAAAAAAGCAGGCTATGCCGTAAACCCAGCTTTCTTGTACAAAGTGGT"
    target = donor[:651] + insert[8:51] + donor[2902:]
    reciprocal = insert[:8] + donor[651:2902] + insert[51:]
    circular = next(p["molecule"]["sequence"]["bases"] for p in products if p["molecule"]["topology"] == "Circular")
    linear = next(p["molecule"]["sequence"]["bases"] for p in products if p["molecule"]["topology"] == "Linear")
    assert len(circular) == len(target) and circular in target + target
    assert linear == reciprocal
    browser.click("Check independent target")
    browser.wait("!!document.querySelector('.recombination-target-check')")
    target_text = browser.evaluate("document.querySelector('.recombination-target-check').textContent")
    assert "1 matches" in target_text and "search complete: true" in target_text, target_text
    assert browser.evaluate("document.querySelectorAll('input[name=recombination-component]:checked').length===0"), "Target matching selected a component implicitly"
    # A user can select the reciprocal component, but it cannot pass as the
    # independently required circular target.
    select_length(browser, 2267)
    verify_selected(browser, "Invalid")
    assert "does not match" in browser.evaluate("document.querySelector('.verification-result').textContent")
    panel(browser); configure(browser); compute(browser)
    select_length(browser, 2553)
    verify_selected(browser, "Incomplete")
    text = browser.evaluate("document.querySelector('.verification-result').textContent")
    assert "matches" in text.lower() and "reference" in text.lower(), text
    html_path = DIRECTORY / "cloning-protocol.html"
    html_path.unlink(missing_ok=True)
    browser.click("Export HTML protocol")
    html = wait_file(html_path).decode()
    assert "Incomplete" in html and "2026-11-24.1" in html and "BP_insert" in html and "BP_donor" in html
    pdf_path = DIRECTORY / "cloning-protocol.pdf"
    pdf_path.unlink(missing_ok=True)
    browser.click("Export PDF protocol")
    pdf = wait_file(pdf_path)
    assert pdf.startswith(b"%PDF-")
    # Independent development validator; the application renders this PDF in Rust.
    pdf_text = subprocess.run(
        ["pdftotext", str(pdf_path), "-"], check=True, capture_output=True, text=True
    ).stdout
    assert "DRAFT - INCOMPLETE" in pdf_text, pdf_text
    assert all(term in pdf_text for term in ["2026-11-24.1", "BP_insert", "BP_donor"]), pdf_text
    # The other product is also a legitimate independently chosen requirement.
    active(browser, "BP_reciprocal"); panel(browser); configure(browser); compute(browser)
    select_length(browser, 2267)
    verify_selected(browser, "Incomplete")
    assert "matches" in browser.evaluate("document.querySelector('.verification-result').textContent").lower()
    # Changing the reaction invalidates both the product evidence and target check.
    panel(browser); configure(browser); compute(browser)
    fill(browser, "Recombination reaction", "LR", "change")
    browser.wait("document.querySelectorAll('input[name=recombination-component]').length===0")
    assert not browser.evaluate("!!document.querySelector('.recombination-target-check')")
    assert not browser.exceptions, browser.exceptions
    assert not browser.console_errors, browser.console_errors
    browser.screenshot(DIRECTORY / "desktop.png")
    print("PASS reciprocal BP UI, source oracle, explicit selection, target check, stale invalidation and incomplete HTML/PDF")
    browser.socket.close()


if __name__ == "__main__":
    main()
