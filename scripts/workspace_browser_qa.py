#!/usr/bin/env python3
"""Development QA against isolated Chrome :9224 and the local build on :8080.

This deliberately clears this test origin's IndexedDB. It exercises real files,
worker validation and browser transactions; it is not an application dependency.
"""
import base64
import io
import json
import time
import zipfile
from pathlib import Path

from browser_qa import Browser


ORIGIN = "http://127.0.0.1:8080"
DIRECTORY = Path("/tmp/collagene2-workspace-qa")


def wait_file(path):
    deadline = time.monotonic() + 20
    while time.monotonic() < deadline:
        if path.exists() and not path.with_suffix(path.suffix + ".crdownload").exists():
            return path.read_bytes()
        time.sleep(.1)
    raise AssertionError(f"Download did not finish: {path}")


def manifest(data):
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        return json.loads(archive.read("manifest.json"))


def tamper_blob(data):
    # Recompute ZIP CRC through an independent writer. Only the content-addressed
    # digest is stale, so passing ZIP validation must not make the archive trusted.
    output = io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(data)) as original, zipfile.ZipFile(output, "w", compression=zipfile.ZIP_STORED) as changed:
        damaged = False
        for entry in original.infolist():
            value = original.read(entry.filename)
            if not damaged and entry.filename.startswith("blobs/") and value:
                value = bytes([value[0] ^ 1]) + value[1:]
                damaged = True
            changed.writestr(entry.filename, value)
    assert damaged
    return output.getvalue()


def slots(browser):
    return browser.evaluate("""new Promise((resolve,reject)=>{
      const request=indexedDB.open('collagene',1);
      request.onerror=()=>reject(Error('Cannot open QA database'));
      request.onsuccess=()=>{
        const db=request.result,tx=db.transaction('workspace'),store=tx.objectStore('workspace'),out={};
        for(const key of ['current','previous']){
          const read=store.get(key);read.onsuccess=()=>{
            if(read.result===undefined){out[key]=null;return;}
            const bytes=new Uint8Array(read.result);let encoded='';
            for(let i=0;i<bytes.length;i+=8192)encoded+=String.fromCharCode(...bytes.subarray(i,i+8192));
            out[key]=btoa(encoded);
          };
        }
        tx.oncomplete=()=>{db.close();resolve(out)};
        tx.onabort=()=>reject(Error('QA read aborted'));
      };
    })""")


def main():
    b = Browser()
    DIRECTORY.mkdir(exist_ok=True)
    for path in DIRECTORY.iterdir():
        if path.is_file():
            path.unlink()
    b.call("Browser.setDownloadBehavior", {"behavior": "allow", "downloadPath": str(DIRECTORY)})
    b.call("Emulation.setDeviceMetricsOverride", {"width": 1600, "height": 1000, "deviceScaleFactor": 1, "mobile": False})
    b.call("Storage.clearDataForOrigin", {"origin": ORIGIN, "storageTypes": "indexeddb"})
    b.call("Page.navigate", {"url": ORIGIN})
    b.wait("document.querySelector('h1')?.textContent==='pCG_demo'")

    def saved():
        b.wait("document.body.innerText.includes('Saved in this browser')")

    def title(name):
        b.wait("document.querySelector('h1')?.textContent===" + json.dumps(name))

    def fill(selector, value):
        b.evaluate("(()=>{const e=document.querySelector(" + json.dumps(selector) + ");e.value=" + json.dumps(value) + ";e.dispatchEvent(new Event('input',{bubbles:true}));})()")

    def close_dialog():
        b.evaluate("document.querySelector('button[aria-label=\"Close dialog\"]').click()")
        b.wait("!document.querySelector('[role=dialog]')")

    def open_file(path, workspace=False):
        selector = '[role=dialog] input[type=file][accept*=cg2workspace]' if workspace else 'input[type=file]'
        b.call("DOM.enable")
        root = b.call("DOM.getDocument")["root"]["nodeId"]
        node = b.call("DOM.querySelector", {"nodeId": root, "selector": selector})["nodeId"]
        assert node, selector
        b.call("DOM.setFileInputFiles", {"nodeId": node, "files": [str(path)]})

    def choose_document(name):
        b.evaluate("Array.from(document.querySelectorAll('.document-entry')).find(e=>e.querySelector('strong').textContent===" + json.dumps(name) + ").click()")
        title(name)

    def rename(name):
        b.click("Document settings")
        fill('input[aria-label="Document name"]', name)
        b.click("Apply document settings")
        title(name)
        b.wait("!document.querySelector('[role=dialog]')")

    def navigate_history():
        b.click("History")
        b.wait("!!document.querySelector('.history-panel')")

    def current_event():
        navigate_history()
        return b.evaluate("Array.from(document.querySelectorAll('.history-panel > details code')).at(-1)?.textContent")

    def events():
        return b.evaluate("Array.from(document.querySelectorAll('.history-card')).map(e=>e.querySelector('details code')?.textContent).filter(e=>e?.startsWith('event-'))")

    def restore_event(event, expected_title):
        navigate_history()
        b.evaluate("(()=>{const card=Array.from(document.querySelectorAll('.history-card')).find(e=>e.querySelector('details code')?.textContent===" + json.dumps(event) + ");Array.from(card.querySelectorAll('button')).find(e=>e.textContent==='Restore revision').click()})()")
        title(expected_title)
        assert current_event() == event

    def undo(expected_title):
        assert b.evaluate("!document.querySelector('button[title=Undo]').disabled")
        b.evaluate("document.querySelector('button[title=Undo]').click()")
        title(expected_title)

    def redo(expected_title):
        assert b.evaluate("!document.querySelector('button[title=Redo]').disabled")
        b.evaluate("document.querySelector('button[title=Redo]').click()")
        title(expected_title)

    saved()
    b.click("Paste sequence")
    fill('textarea[aria-label="Sequence input"]', '>left\nAAAACCCC\n>right\nCCCCGGGG\n>target\nAAAACCCCGGGG\n')
    b.click("Open sequence")
    title("left")
    choose_document("target")
    root_event = current_event()
    assert root_event.startswith("event-")
    rename("target_A")
    event_a = current_event()
    rename("target_B")
    event_b = current_event()
    undo("target_A")
    saved()
    b.call("Page.reload")
    title("target_A")
    assert current_event() == event_a
    redo("target_B")
    assert current_event() == event_b
    undo("target_A")
    rename("target_C")
    event_c = current_event()
    assert b.evaluate("document.querySelector('button[title=Redo]').disabled")
    assert set(events()) == {root_event, event_a, event_b, event_c}
    assert b.evaluate("""(()=>{
      const panel=document.querySelector('.history-panel');
      const next=Array.from(panel.querySelectorAll('button')).filter(e=>e.textContent.trim()==='Next 100 events');
      return next.length===1 && next[0].disabled
        && !panel.textContent.includes('events.with(Vec::len)')
        && !panel.textContent.includes('on : click');
    })()"""), "History pagination must show a disabled Next 100 events button without leaked Rust source"
    restore_event(event_b, "target_B")
    undo("target_C")
    redo("target_B")
    restore_event(root_event, "target")
    saved()
    b.call("Page.reload")
    title("target")
    navigate_history()
    assert set(events()) == {root_event, event_a, event_b, event_c}

    # Build a candidate in the real worker, then import its native exchange file.
    b.click("Plan")
    b.evaluate("""(()=>{
      for(const [label,name] of [['First source document','left'],['Second source document','right']]){
        const s=document.querySelector('select[aria-label="'+label+'"]');
        s.value=Array.from(s.options).find(o=>o.textContent===name).value;
        s.dispatchEvent(new Event('change',{bubbles:true}));
      }
      const i=Array.from(document.querySelectorAll('.inspector label')).find(e=>e.textContent.includes('Junction overlap')).querySelector('input');
      i.value='4';i.dispatchEvent(new Event('input',{bubbles:true}));
    })()""")
    b.click("Reconstruct and verify")
    b.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    b.click("Save candidate")
    candidate_path = DIRECTORY / "candidate.cg2plan"
    wait_file(candidate_path)
    open_file(candidate_path)
    b.wait("!!document.querySelector('.candidate-review')")
    assert not b.evaluate("!!document.querySelector('.verification-result')")
    saved()
    b.click("Workspace")
    b.click("Save workspace file")
    archive_path = DIRECTORY / "collagene-workspace.cg2workspace"
    archive = wait_file(archive_path)
    expected = manifest(archive)
    assert expected["active_candidate"]
    assert len(expected["candidates"]) == 1
    assert {e["id"] for e in expected["history"] if e["document_id"] == expected["active_document"]} == {root_event, event_a, event_b, event_c}
    expected_document = next(d for d in expected["documents"] if d["id"] == expected["active_document"])
    assert expected_document["head_event"] == root_event
    assert expected_document["undo_events"]
    close_dialog()

    # Replacement makes a recoverable previous workspace, including a deliberately
    # closed candidate: retaining it must not resurrect its active marker.
    choose_document("left")
    rename("replacement_workspace")
    b.click("Plan")
    b.click("Discard candidate")
    saved()
    previous_state = manifest(base64.b64decode(slots(b)["current"]))
    assert previous_state["active_document"] != expected["active_document"]
    assert previous_state["active_candidate"] is None
    assert len(previous_state["candidates"]) == 1
    b.click("Workspace")
    open_file(archive_path, workspace=True)
    title("target")
    b.wait("!document.querySelector('[role=dialog]')")
    b.click("Plan")
    b.wait("!!document.querySelector('.candidate-review')")
    assert not b.evaluate("!!document.querySelector('.verification-result')"), "An archive must not import a trusted verification result"
    saved()
    stored = slots(b)
    restored = manifest(base64.b64decode(stored["current"]))
    assert restored["documents"] == expected["documents"]
    assert restored["history"] == expected["history"]
    assert restored["active_document"] == expected["active_document"]
    assert restored["active_candidate"] == expected["active_candidate"]
    assert manifest(base64.b64decode(stored["previous"]))["documents"] == previous_state["documents"]

    # A ZIP with good CRCs and stale content hashes changes neither storage slot.
    corrupt_path = DIRECTORY / "tampered.cg2workspace"
    corrupt_path.write_bytes(tamper_blob(archive))
    before_names = b.evaluate("Array.from(document.querySelectorAll('.document-entry strong')).map(e=>e.textContent)")
    b.click("Workspace")
    open_file(corrupt_path, workspace=True)
    b.wait("document.querySelector('[role=dialog]')?.textContent.includes('Workspace was not replaced')")
    assert slots(b) == stored, "Failed archive import changed IndexedDB"
    assert b.evaluate("document.querySelector('h1').textContent") == "target"
    assert b.evaluate("Array.from(document.querySelectorAll('.document-entry strong')).map(e=>e.textContent)") == before_names
    close_dialog()
    b.click("Plan")
    assert b.evaluate("!!document.querySelector('.candidate-review')")
    b.click("Workspace")
    b.click("Restore previous workspace")
    title("replacement_workspace")
    b.wait("!document.querySelector('[role=dialog]')")
    b.click("Plan")
    assert not b.evaluate("!!document.querySelector('.candidate-review')"), "Closed candidate was resurrected"
    saved()
    recovered = manifest(base64.b64decode(slots(b)["current"]))
    assert recovered["documents"] == previous_state["documents"]
    assert recovered["active_candidate"] is None
    navigate_history()
    b.wait("document.querySelector('.history-panel')?.textContent.includes('Open candidate for review')")

    # Corrupt current autosave must remain intact at boot; explicit recovery can
    # still restore the valid previous slot through the same UI validation path.
    b.evaluate("""new Promise((resolve,reject)=>{const r=indexedDB.open('collagene',1);r.onsuccess=()=>{const d=r.result,t=d.transaction('workspace','readwrite');t.objectStore('workspace').put(new Uint8Array([80,75,0,1]),'current');t.oncomplete=()=>{d.close();resolve(true)};t.onabort=()=>reject(Error('QA corruption write failed'))}})""")
    before_boot = slots(b)
    b.call("Page.reload")
    b.wait("document.body?.innerText.includes('Recovery needed')")
    assert slots(b) == before_boot
    b.click("Workspace")
    b.click("Restore previous workspace")
    title("target")
    b.wait("!document.querySelector('[role=dialog]')")
    b.click("Plan")
    b.wait("!!document.querySelector('.candidate-review')")
    assert not b.evaluate("!!document.querySelector('.verification-result')")
    saved()
    b.call("Page.reload")
    title("target")
    navigate_history()
    assert set(events()) == {root_event, event_a, event_b, event_c}
    b.screenshot("/tmp/collagene2-workspace-history.png")
    assert not b.exceptions, json.dumps(b.exceptions)
    assert not b.console_errors, json.dumps(b.console_errors)
    print("Workspace browser QA passed: reload undo/redo, retained branches, archive save/open, active and closed candidate preservation, tampered archive atomic rejection, previous-slot recovery and corrupt-autosave recovery.")


if __name__ == "__main__":
    main()
