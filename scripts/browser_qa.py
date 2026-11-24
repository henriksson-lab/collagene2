"""Development smoke test against an already-running isolated Chrome and local build.

Run Chrome with CDP on 127.0.0.1:9224 and serve dist on 127.0.0.1:8080.
This helper is not part of the application's runtime.
"""
import base64
import json
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import websocket

QA_ORIGIN = os.environ.get("COLLAGENE_QA_ORIGIN", 'http://127.0.0.1:8080').rstrip("/")
QA_DIST = Path(os.environ.get("COLLAGENE_QA_DIST", str(Path(__file__).resolve().parents[1] / "dist")))


def mobile_enabled():
    """Mobile acceptance is deferred; opt in explicitly after desktop review."""
    return os.environ.get("COLLAGENE_QA_MOBILE") == "1"


class Browser:
    def __init__(self):
        targets = json.load(urllib.request.urlopen("http://127.0.0.1:9224/json/list"))
        target = next(t for t in targets if t["type"] == "page")
        self.socket = websocket.create_connection(
            target["webSocketDebuggerUrl"], origin="http://127.0.0.1:9224", timeout=30
        )
        self.serial = 0
        self.exceptions = []
        self.console_errors = []
        self.call("Runtime.enable")
        self.call("Page.enable")
        self.call("Log.enable")

    def call(self, method, params=None):
        self.serial += 1
        serial = self.serial
        self.socket.send(json.dumps({"id": serial, "method": method, "params": params or {}}))
        while True:
            result = json.loads(self.socket.recv())
            if result.get("method") == "Fetch.requestPaused" and (getattr(self,"fetch_fixture",None) or getattr(self,"fetch_routes",None)):
                request=result["params"];url=request["request"]["url"]
                fixture=next((route for route in getattr(self,"fetch_routes",[]) if url.startswith(route["prefix"])),getattr(self,"fetch_fixture",None))
                self.serial+=1
                if fixture and url.startswith(fixture["prefix"]):
                    if fixture.get('defer') and request['request']['method']!='OPTIONS':
                        self.deferred_requests=getattr(self,'deferred_requests',[])+[request['requestId']]
                        continue
                    self.fulfilled_urls=getattr(self,"fulfilled_urls",[])+[url]
                    if request['request']['method']!='OPTIONS' and 'authorization' in fixture:
                        actual=next((v for k,v in request['request'].get('headers',{}).items() if k.lower()=='authorization'),None)
                        if actual!=fixture['authorization']:
                            self.request_violations=getattr(self,'request_violations',[])+['Unexpected authorization routing for mocked source']
                    params={"requestId":request["requestId"],"responseCode":fixture.get('status',200),"responseHeaders":[{"name":"Content-Type","value":fixture.get('mime','text/plain')},{"name":"Access-Control-Allow-Origin","value":"*"},{"name":"Access-Control-Allow-Headers","value":"Authorization"},{"name":"Access-Control-Allow-Methods","value":"GET, OPTIONS"}],"body":base64.b64encode(b'' if request['request']['method']=='OPTIONS' else fixture["body"]).decode("ascii")}
                    method="Fetch.fulfillRequest"
                else:
                    method="Fetch.failRequest";params={"requestId":request["requestId"],"errorReason":"Aborted"}
                self.socket.send(json.dumps({"id":self.serial,"method":method,"params":params}))
            if result.get("method") == "Runtime.exceptionThrown":
                self.exceptions.append(result["params"])
            if result.get("method") == "Log.entryAdded":
                entry=result["params"]["entry"]
                if entry.get("level")=="error" or "ResizeObserver" in entry.get("text",""):
                    self.console_errors.append(entry)
            if result.get("id") == serial:
                if "error" in result:
                    raise RuntimeError(result["error"])
                return result.get("result", {})

    def evaluate(self, expression):
        result = self.call("Runtime.evaluate", {"expression": expression, "returnByValue": True, "awaitPromise": True})
        if "exceptionDetails" in result:
            raise RuntimeError(result["exceptionDetails"])
        return result.get("result", {}).get("value")

    def wait(self, expression, timeout=30):
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if self.evaluate(f"Boolean({expression})"):
                return
            time.sleep(0.1)
        raise AssertionError(f"Timed out: {expression}\n{self.evaluate('document.body.innerText')}\nExceptions: {json.dumps(self.exceptions)}\nConsole: {json.dumps(self.console_errors[-10:])}")

    def pointer_click(self, expression):
        """Require a visible, unobscured pointer target; never bypass closed details."""
        # Event lookup can schedule focus on the next animation frame. Settle that
        # focus before scrolling; otherwise it can scroll the inspector back down
        # after we bring navigation into view. Keep actual visibility/hit checks.
        self.evaluate("new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))")
        self.evaluate(f"({expression}).scrollIntoView({{block:'nearest',inline:'nearest',behavior:'instant'}})")
        self.evaluate("new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))")
        point = self.evaluate("(()=>{const e=(" + expression + ");if(!e)return null;"
            "const r=e.getBoundingClientRect(),s=getComputedStyle(e),x=r.left+r.width/2,y=r.top+r.height/2;"
            "const hit=document.elementFromPoint(x,y);return {x,y,width:r.width,height:r.height,"
            "visible:s.visibility==='visible'&&s.display!=='none'&&r.width>0&&r.height>0,"
            "inViewport:x>=0&&x<innerWidth&&y>=0&&y<innerHeight,"
            "unobscured:!!hit&&(hit===e||e.contains(hit)),disabled:!!e.disabled};})()")
        assert point and point['visible'] and point['inViewport'] and point['unobscured'] and not point['disabled'], (expression, point)
        self.call('Input.dispatchMouseEvent', {'type':'mouseMoved','x':point['x'],'y':point['y']})
        self.call('Input.dispatchMouseEvent', {'type':'mousePressed','x':point['x'],'y':point['y'],'button':'left','clickCount':1})
        self.call('Input.dispatchMouseEvent', {'type':'mouseReleased','x':point['x'],'y':point['y'],'button':'left','clickCount':1})

    def click(self, label):
        text = json.dumps(label)
        # Prefer registry navigation over a duplicate toolbar shortcut, so every
        # secondary tool is proven discoverable through the visible More tools menu.
        tool = "Array.from(document.querySelectorAll('.more-tools [data-tool-id]')).find(b=>b.textContent.trim()===" + text + ")"
        if self.evaluate(f"!!({tool})"):
            tool_id = self.evaluate(f"({tool}).dataset.toolId")
            assert tool_id, label
            menu = f"({tool}).closest('details')"
            if not self.evaluate(f"({menu}).open"):
                self.pointer_click(f"({menu}).querySelector('summary[aria-label=\"More tools\"]')")
                self.wait(f"({menu}).open")
            self.pointer_click(tool)
            self.wait(f"!({menu}).open")
            return
        primary = "Array.from(document.querySelectorAll('.inspector-tabs [data-panel-id]')).find(b=>b.textContent.trim()===" + text + ")"
        if self.evaluate(f"!!({primary})"):
            self.pointer_click(primary)
            return
        expression = "Array.from(document.querySelectorAll('button')).find(b=>b.textContent.trim()===" + text + ")"
        assert self.evaluate(f"!!({expression})"), label
        self.evaluate(f"({expression}).click()")

    def screenshot(self, path):
        image = self.call("Page.captureScreenshot", {"format": "png", "captureBeyondViewport": False})
        Path(path).write_bytes(base64.b64decode(image["data"]))


def main():
    with urllib.request.urlopen(QA_ORIGIN,timeout=3) as response:
        assert response.status==200, "Local preview is not ready"
    browser = Browser()
    downloads=Path("/tmp/collagene2-qa-downloads")
    downloads.mkdir(exist_ok=True)
    for name in ["target.cg2","cloning-protocol.html","cloning-protocol.pdf","planning-workspace.cg2design","candidate.cg2plan"]:
        (downloads/name).unlink(missing_ok=True)
    browser.call("Browser.setDownloadBehavior", {"behavior":"allow","downloadPath":str(downloads)})
    browser.call("Emulation.setDeviceMetricsOverride", {"width": 1600, "height": 1000, "deviceScaleFactor": 1, "mobile": False})
    browser.call("Storage.clearDataForOrigin", {"origin":QA_ORIGIN,"storageTypes":"indexeddb"})
    browser.call("Page.navigate", {"url": QA_ORIGIN})
    browser.wait("document.querySelector('h1')?.textContent==='pCG_demo'")
    browser.wait("document.body.innerText.includes('Saved in this browser')")
    browser.wait("document.querySelectorAll('canvas').length===3")
    assert not browser.evaluate("document.body.innerText.includes('WebGL2 is unavailable')")
    browser.evaluate("new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))")
    assert browser.evaluate("document.querySelector('.map-labels').getBoundingClientRect().height<1000"), "Map grew outside its viewport"
    assert browser.evaluate("(()=>{const c=document.querySelector('.map-labels');return c.getContext('2d').getImageData(0,0,c.width,c.height).data.some(v=>v!==0)})()"), "Map labels were not painted"
    browser.screenshot("/tmp/collagene2-desktop.png")

    browser.evaluate("Array.from(document.querySelectorAll('.feature-entry')).find(b=>b.textContent.includes('Reporter example')).click()")
    browser.wait("document.querySelector('.statusbar').textContent.includes('451–1250')")
    browser.click("Edit selection")
    original=browser.evaluate("document.querySelector('textarea').value")
    browser.evaluate("(()=>{const t=document.querySelector('textarea');t.value='?';t.dispatchEvent(new Event('input',{bubbles:true}));})()")
    browser.click("Apply edit")
    browser.wait("document.querySelector('[role=dialog] [role=alert]')?.textContent.includes('Invalid')")
    replacement=('A' if original[0]!='A' else 'C')+original[1:]
    browser.evaluate("(()=>{const t=document.querySelector('textarea');t.value="+json.dumps(replacement)+";t.dispatchEvent(new Event('input',{bubbles:true}));})()")
    browser.click("Apply edit")
    browser.wait("document.querySelector('.feature-list')?.textContent.includes('review required')")
    browser.evaluate("document.querySelector('button[title=Undo]').click()")
    browser.wait("!document.querySelector('.feature-list')?.textContent.includes('review required')")

    browser.click("Digest")
    assert browser.evaluate('document.querySelector(\'[aria-label="Digest enzyme EcoRI"]\').checked')
    browser.click("Run virtual digest")
    browser.wait("!!document.querySelector('.digest-results[data-scan-complete=true][data-products-determined=true]')")
    assert browser.evaluate("document.querySelector('.fragment-sizes').textContent.includes('nt')")
    assert browser.evaluate("Array.from(document.querySelector('.fragment-sizes').textContent.matchAll(/(\\d+) nt/g)).reduce((sum,m)=>sum+Number(m[1]),0)===4800")

    browser.evaluate("(()=>{const i=document.querySelector('.search input');i.value='GAATTC';i.dispatchEvent(new Event('input',{bubbles:true}));})()")
    browser.click("Find")
    browser.wait("document.querySelector('.notice')?.textContent.includes('match(es)')")

    browser.click("Paste sequence")
    browser.evaluate("(()=>{const t=document.querySelector('textarea');t.value='?';t.dispatchEvent(new Event('input',{bubbles:true}));})()")
    browser.click("Open sequence")
    browser.wait("document.querySelector('[role=dialog] [role=alert]')?.textContent.includes('Invalid')")
    sequence = ">left\nAAAACCCC\n>right\nCCCCGGGG\n>target\nAAAACCCCGGGG\n>pcr_template\nACGACCCCGGGGTAAC\n"
    browser.evaluate("(()=>{const t=document.querySelector('textarea');t.value=" + json.dumps(sequence) + ";t.dispatchEvent(new Event('input',{bubbles:true}));})()")
    browser.click("Open sequence")
    browser.wait("document.querySelector('h1')?.textContent==='left'")
    browser.evaluate("Array.from(document.querySelectorAll('.document-entry')).find(b=>b.textContent.includes('target')).click()")
    browser.wait("document.querySelector('h1')?.textContent==='target'")
    browser.click("Plan")
    browser.evaluate("(()=>{const selects=document.querySelectorAll('.inspector select');for(const [i,value] of [[0,'1'],[1,'2']]){selects[i].value=value;selects[i].dispatchEvent(new Event('change',{bubbles:true}));}const input=document.querySelector('.inspector input');input.value='4';input.dispatchEvent(new Event('input',{bubbles:true}));})()")
    browser.click("Reconstruct and verify")
    browser.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    browser.click("Export HTML protocol")
    browser.click("Export PDF protocol")
    browser.click("Save .cg2")
    for name in ["cloning-protocol.html","cloning-protocol.pdf","target.cg2"]:
        deadline=time.monotonic()+10
        while not (downloads/name).exists() and time.monotonic()<deadline:
            time.sleep(0.1)
        assert (downloads/name).exists(), name
    assert (downloads/"target.cg2").read_bytes().startswith(b"CG2\0")
    assert (downloads/"cloning-protocol.pdf").read_bytes().startswith(b"%PDF-")
    extracted=subprocess.run(["pdftotext",str(downloads/"cloning-protocol.pdf"),"-"],capture_output=True,text=True,check=True).stdout
    assert "Consistent under assumptions" in extracted and "target" in extracted
    protocol_html=(downloads/"cloning-protocol.html").read_text()
    assert "ConsistentUnderAssumptions" in protocol_html and "Reconstructed product matches" in protocol_html

    # Deliberately wrong overlap must fail, never reuse the earlier successful report.
    browser.evaluate("(()=>{const input=document.querySelector('.inspector input');input.value='3';input.dispatchEvent(new Event('input',{bubbles:true}));})()")
    browser.wait("!document.querySelector('.verification-result')")
    browser.click("Reconstruct and verify")
    browser.wait("document.querySelector('.verification-result')?.textContent.includes('Invalid')")

    # A saved snapshot survives page reload and reopens the same target.
    browser.wait("document.body.innerText.includes('Saved in this browser')")

    # Verify native file opening through the actual browser File API.
    browser.call("DOM.enable")
    root=browser.call("DOM.getDocument")["root"]["nodeId"]
    file_input=browser.call("DOM.querySelector", {"nodeId":root,"selector":"input[type=file]"})["nodeId"]
    before=browser.evaluate("document.querySelectorAll('.document-entry').length")
    browser.call("DOM.setFileInputFiles", {"nodeId":file_input,"files":[str(downloads/"target.cg2")]})
    browser.wait(f"document.querySelectorAll('.document-entry').length==={before+1}")
    browser.wait("document.querySelector('h1')?.textContent==='target'")

    # Browser snapshot -> actual stdio MCP server -> binary candidate -> WASM verification.
    browser.click("Plan")
    browser.click("Export for AI")
    workspace_file=downloads/"planning-workspace.cg2design"
    deadline=time.monotonic()+10
    while not workspace_file.exists() and time.monotonic()<deadline:
        time.sleep(0.1)
    assert workspace_file.read_bytes().startswith(b"CGD\0")
    root_path=Path(__file__).resolve().parent.parent
    subprocess.run([sys.executable,str(root_path/"scripts/mcp_smoke.py"),"--workspace",str(workspace_file),"--candidate",str(downloads/"candidate.cg2plan")],cwd=root_path,check=True,timeout=60)
    browser.call("DOM.setFileInputFiles", {"nodeId":file_input,"files":[str(downloads/"candidate.cg2plan")]})
    browser.wait("document.querySelector('.candidate-review')?.textContent.includes('Stdio assembly')")
    assert not browser.evaluate("!!document.querySelector('.verification-result')"), "Imported candidates must not carry trusted results"
    assert browser.evaluate("document.querySelector('h1').textContent==='target'"), "Candidate import replaced the target"
    browser.click("Verify candidate")
    browser.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    (downloads/"cloning-protocol.html").unlink()
    browser.click("Export HTML protocol")
    deadline=time.monotonic()+10
    while not (downloads/"cloning-protocol.html").exists() and time.monotonic()<deadline:
        time.sleep(0.1)
    protocol_html=(downloads/"cloning-protocol.html").read_text()
    assert "Stdio assembly" in protocol_html and "ConsistentUnderAssumptions" in protocol_html

    def choose_document(name):
        browser.evaluate("Array.from(document.querySelectorAll('.document-entry')).find(b=>b.textContent.includes("+json.dumps(name)+")).click()")
        browser.wait("document.querySelector('h1')?.textContent==="+json.dumps(name))

    def edit_first_base():
        browser.evaluate("(()=>{const c=document.querySelector('.sequence-canvas');const r=c.getBoundingClientRect();c.dispatchEvent(new MouseEvent('click',{bubbles:true,clientX:r.left+68,clientY:r.top+30}));})()")
        browser.click("Edit selection")
        browser.evaluate("(()=>{const t=document.querySelector('textarea');t.value='T';t.dispatchEvent(new Event('input',{bubbles:true}));})()")
        browser.click("Apply edit")
        browser.wait("!document.querySelector('[role=dialog]')")

    choose_document("left")
    edit_first_base()
    choose_document("target")
    browser.click("Verify candidate")
    browser.wait("document.querySelector('.notice')?.textContent.includes('missing or has changed')")
    assert not browser.evaluate("!!document.querySelector('.verification-result')")
    choose_document("left")
    browser.evaluate("document.querySelector('button[title=Undo]').click()")
    choose_document("target")
    edit_first_base()
    browser.click("Verify candidate")
    browser.wait("document.querySelector('.notice')?.textContent.includes('active target differs')")
    browser.evaluate("document.querySelector('button[title=Undo]').click()")
    browser.click("Verify candidate")
    browser.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    browser.wait("document.body.innerText.includes('Saved in this browser')")
    browser.call("Page.reload")
    browser.wait("document.querySelector('h1')?.textContent==='target'")
    browser.click("Plan")
    browser.wait("document.querySelector('.candidate-review')?.textContent.includes('Stdio assembly')")
    assert not browser.evaluate("!!document.querySelector('.verification-result')"), "Recovery must reverify candidates"
    browser.click("Verify candidate")
    browser.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    browser.screenshot("/tmp/collagene2-candidate.png")
    browser.click("Discard candidate")

    # Independent restriction/ligation fixture, including explicitly retained nicks.
    root=browser.call("DOM.getDocument")["root"]["nodeId"]
    file_input=browser.call("DOM.querySelector", {"nodeId":root,"selector":"input[type=file]"})["nodeId"]
    browser.call("DOM.setFileInputFiles", {"nodeId":file_input,"files":[str(root_path/"fixtures/restriction-demo.gb")]})
    browser.wait("document.querySelector('h1')?.textContent==='vector_fixture'")
    browser.click("Digest")
    assert browser.evaluate('document.querySelector(\'[aria-label="Digest enzyme EcoRI"]\').checked')
    browser.click("Run virtual digest")
    browser.wait("!!document.querySelector('.digest-results[data-scan-complete=true][data-products-determined=true]')")
    assert browser.evaluate("document.querySelectorAll('.digest-fragment').length===1")
    assert browser.evaluate("document.querySelector('.fragment-sizes').textContent==='14 nt'")
    assert browser.evaluate("document.querySelector('.digest-fragment').textContent.includes('14 nt top / 14 nt bottom')")
    assert browser.evaluate("document.querySelector('.digest-fragment').textContent.includes('original reference start 6')")
    assert browser.evaluate("document.querySelector('.digest-fragment').textContent.includes('Left: 5′ overhang AATT · right: 5′ overhang AATT')")
    choose_document("ligation_target")
    browser.click("Plan")
    browser.evaluate("(()=>{for(const [label,name] of [['First source document','vector_fixture'],['Second source document','insert_fixture']]){const s=document.querySelector('select[aria-label=\"'+label+'\"]');s.value=Array.from(s.options).find(o=>o.textContent===name).value;s.dispatchEvent(new Event('change',{bubbles:true}));}const method=document.querySelector('select[aria-label=\"Cloning method\"]');method.value='restriction';method.dispatchEvent(new Event('change',{bubbles:true}));})()")
    browser.wait("!!document.querySelector('input[aria-label=\"First enzymes\"]')")
    browser.click("Reconstruct and verify")
    browser.wait("document.querySelector('.verification-result')?.textContent.includes('ConsistentUnderAssumptions')")
    assert browser.evaluate("document.querySelector('.verification-result').textContent.includes('0 retained nicks')")
    browser.evaluate("(()=>{const t=document.querySelector('select[aria-label=\"First end treatment\"]');t.value='dephosphorylate';t.dispatchEvent(new Event('change',{bubbles:true}));})()")
    browser.wait("!document.querySelector('.verification-result')")
    browser.click("Reconstruct and verify")
    browser.wait("document.querySelector('.verification-result')?.textContent.includes('Invalid')")
    assert browser.evaluate("document.querySelector('.verification-result').textContent.includes('requires sealed junctions')")
    browser.evaluate("(()=>{const c=Array.from(document.querySelectorAll('label')).find(l=>l.textContent.includes('Allow explicitly reported nicks')).querySelector('input');c.checked=true;c.dispatchEvent(new Event('change',{bubbles:true}));})()")
    browser.click("Reconstruct and verify")
    browser.wait("document.querySelector('.verification-result')?.textContent.includes('2 retained nicks')")
    assert browser.evaluate("document.querySelector('.verification-result').textContent.includes('Covalently closed in this model: false')")
    (downloads/"cloning-protocol.html").unlink()
    browser.click("Export HTML protocol")
    deadline=time.monotonic()+10
    while not (downloads/"cloning-protocol.html").exists() and time.monotonic()<deadline:
        time.sleep(0.1)
    protocol_html=(downloads/"cloning-protocol.html").read_text()
    assert "2 retained nicks" in protocol_html and "Complete restriction digest" in protocol_html and "Ordered ligation" in protocol_html
    choose_document("target")

    # Print the exported protocol as the browser PDF path, then return to the app.
    frame=browser.call("Page.getFrameTree")["frameTree"]["frame"]["id"]
    browser.wait("document.body?.innerText.includes('Saved in this browser')")
    browser.call("Page.setDocumentContent", {"frameId":frame,"html":protocol_html})
    pdf=base64.b64decode(browser.call("Page.printToPDF", {"printBackground":True})["data"])
    assert pdf.startswith(b"%PDF-") and len(pdf)>1000
    Path("/tmp/collagene2-protocol.pdf").write_bytes(pdf)
    browser.call("Page.navigate", {"url":QA_ORIGIN})
    browser.wait("document.querySelector('h1')?.textContent==='target'")
    # Run the same manually specified tailed-PCR fixture as the native scientific test.
    browser.evaluate("Array.from(document.querySelectorAll('.document-entry')).find(b=>b.textContent.includes('pcr_template')).click()")
    browser.wait("document.querySelector('h1')?.textContent==='pcr_template'")
    browser.click("PCR")
    browser.evaluate("(()=>{const inputs=document.querySelectorAll('.inspector input');['ACGA','CGT','GTTA','GAC'].forEach((value,i)=>{inputs[i].value=value;inputs[i].dispatchEvent(new Event('input',{bubbles:true}));});})()")
    browser.click("Compute products")
    browser.wait("document.querySelector('h1')?.textContent==='pcr_template_PCR_1'")
    assert browser.evaluate("document.querySelector('.document-meta').textContent.includes('22 bp')")
    browser.evaluate("(()=>{const c=document.querySelector('.sequence-canvas');const r=c.getBoundingClientRect();for(const [i,shift] of [[0,false],[21,true]]){c.dispatchEvent(new MouseEvent('click',{bubbles:true,clientX:r.left+64+i*8.5+4,clientY:r.top+30,shiftKey:shift}));}})()")
    browser.wait("document.querySelector('.sr-only')?.textContent.includes('CGTACGACCCCGGGGTAACGTC')")
    browser.evaluate("Array.from(document.querySelectorAll('.document-entry')).find(b=>b.textContent.includes('target')).click()")
    browser.wait("document.querySelector('h1')?.textContent==='target'")
    browser.wait("document.body?.innerText.includes('Saved in this browser')")
    browser.call("Page.reload")
    browser.wait("document.querySelector('h1')?.textContent==='target'")
    browser.wait("document.body.innerText.includes('Saved in this browser')")

    if mobile_enabled():
        browser.call("Emulation.setDeviceMetricsOverride", {"width": 390, "height": 844, "deviceScaleFactor": 1, "mobile": True})
    browser.call("Page.reload")
    browser.wait("document.querySelector('h1')?.textContent==='target'")
    if mobile_enabled():
        assert browser.evaluate("document.documentElement.scrollWidth<=window.innerWidth"), "Mobile layout overflows"
        browser.screenshot("/tmp/collagene2-mobile.png")
    assert not browser.exceptions, json.dumps(browser.exceptions)
    assert not browser.console_errors, json.dumps(browser.console_errors[-10:])
    print("Browser smoke checks passed: actual map paint, IndexedDB recovery, selection/edit/undo and validation errors, worker digest/search/PCR, correct and incorrect assembly verification, native file reopening, browser/MCP binary exchange with stale-input/target rejection and candidate recovery, duplex digest and ligation with explicit nicks, HTML/embedded and browser PDF export, and mobile layout.".replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__ == "__main__":
    import sys
    if "--close" in sys.argv:
        browser=Browser()
        browser.socket.send(json.dumps({"id":999,"method":"Browser.close"}))
        print("Closed isolated test browser.")
    elif "--stop-preview" in sys.argv:
        import os
        import signal
        root=Path(__file__).resolve().parent.parent
        expected_prefix=[b"-m",b"http.server",b"8080",b"--bind",b"127.0.0.1",b"--directory"]
        allowed_directories={b"dist",os.fsencode(root/"dist")}
        for process in Path("/proc").iterdir():
            if not process.name.isdigit():
                continue
            try:
                arguments=(process/"cmdline").read_bytes().rstrip(b"\0").split(b"\0")
                if process.stat().st_uid==os.getuid() and arguments[1:-1]==expected_prefix and arguments[-1] in allowed_directories and (process/"cwd").resolve()==root:
                    os.kill(int(process.name),signal.SIGTERM)
                    print("Stopped this checkout's temporary preview server.")
            except (OSError,ProcessLookupError):
                pass
    elif "--inspect" in sys.argv:
        browser=Browser()
        browser.call("Emulation.setDeviceMetricsOverride", {"width":1600,"height":1000,"deviceScaleFactor":1,"mobile":False})
        browser.call("Page.reload")
        browser.wait("document.querySelector('.map-stage canvas')")
        browser.evaluate("new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))")
        browser.screenshot("/tmp/collagene2-diagnostic.png")
        print(browser.call("Browser.getVersion"))
        print(browser.evaluate("JSON.stringify({height:innerHeight,width:innerWidth,dvh:CSS.supports('height','100dvh'),styles:Array.from(document.styleSheets).map(s=>s.href)})"))
        print(browser.evaluate("JSON.stringify(Array.from(document.querySelectorAll('.application,.body-layout,main,.editor-grid')).map(e=>({class:e.className,height:getComputedStyle(e).height,min:getComputedStyle(e).minHeight,display:getComputedStyle(e).display,flex:getComputedStyle(e).flex})))"))
        print(browser.evaluate("JSON.stringify(Array.from(document.querySelectorAll('.editor-grid,.sequence-panel,.sequence-scroll,.map-panel,.map-stage,canvas')).map(e=>({class:e.className,rect:e.getBoundingClientRect().toJSON(),width:e.width,height:e.height})))"))
    else:
        main()
