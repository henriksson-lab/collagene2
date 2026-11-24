#!/usr/bin/env python3
"""Development QA for real event-DAG rendering/navigation and cross-document forks.

Requires isolated Chrome :9224 and this checkout's release preview :8080.
Clears IndexedDB for that test origin. No application runtime dependency.
"""
import json
from browser_qa import Browser, mobile_enabled


def main():
    b=Browser()
    b.call('Emulation.setDeviceMetricsOverride',{'width':1600,'height':1000,'deviceScaleFactor':1,'mobile':False})
    b.call('Storage.clearDataForOrigin',{'origin':'http://127.0.0.1:8080','storageTypes':'indexeddb'})
    b.call('Page.navigate',{'url':'http://127.0.0.1:8080'})
    b.wait("document.querySelector('h1')?.textContent==='pCG_demo'")
    b.wait("document.body.innerText.includes('Saved in this browser')")

    def fill(selector,value):
        b.evaluate("(()=>{const e=document.querySelector("+json.dumps(selector)+");e.value="+json.dumps(value)+";e.dispatchEvent(new Event('input',{bubbles:true}));})()")

    def history():
        b.click('History')
        b.wait("!!document.querySelector('.history-graph-event')")

    def current():
        history()
        return b.evaluate("document.querySelector('.history-graph-event[data-current=true]').dataset.eventId")

    def title(name):
        b.wait("document.querySelector('h1')?.textContent==="+json.dumps(name))

    def rename(name):
        b.click('Document settings')
        fill('input[aria-label="Document name"]',name)
        b.click('Apply document settings')
        title(name)
        b.wait("!document.querySelector('[role=dialog]')")

    def node(event):
        return 'document.querySelector('+json.dumps('.history-graph-event[data-event-id="'+event+'"]')+')'

    def inspect(event):
        b.evaluate(node(event)+'.click()')
        b.wait("document.querySelector('.history-graph-inspection')?.dataset.eventId==="+json.dumps(event))

    def selected():
        return b.evaluate("document.querySelector('.history-graph-inspection').dataset.eventId")

    def key(name):
        b.call('Input.dispatchKeyEvent',{'type':'keyDown','key':name,'code':name})
        b.call('Input.dispatchKeyEvent',{'type':'keyUp','key':name,'code':name})

    b.click('Paste sequence')
    fill('textarea[aria-label="Sequence input"]','>history_A\nACGTACGT\n')
    b.click('Open sequence')
    title('history_A')
    root=current()
    initial_revision=b.evaluate(node(root)+'.dataset.revision')
    rename('history_B');middle=current()
    rename('history_A');returned=current()
    assert returned!=root
    assert b.evaluate(node(returned)+'.dataset.revision')==initial_revision
    assert b.evaluate(node(root)+".textContent.includes('Same revision · 2 events')")
    assert b.evaluate(node(returned)+".textContent.includes('Same revision · 2 events')")
    assert b.evaluate("document.querySelectorAll('.history-graph-event').length===3")
    # Undo/Redo use exact event stacks despite equal-content revisions.
    b.evaluate("document.querySelector('button[title=Undo]').click()")
    title('history_B');assert current()==middle
    assert b.evaluate(node(returned)+'.dataset.redoDistance')=='1'
    assert b.evaluate(node(root)+'.dataset.undoDistance')=='1'
    rename('history_C');branch=current()
    assert b.evaluate("document.querySelector('button[title=Redo]').disabled")
    assert b.evaluate("document.querySelectorAll('.history-graph-event').length===4")
    assert b.evaluate(node(middle)+".textContent.includes('Branch point')")
    for child in [returned,branch]:
        assert b.evaluate("Array.from(document.querySelectorAll('.history-graph-edge')).some(e=>e.dataset.child==="+json.dumps(child)+"&&e.dataset.parent==="+json.dumps(middle)+")")
    assert not b.evaluate("Array.from(document.querySelectorAll('.history-graph-edge')).some(e=>e.dataset.child==="+json.dumps(branch)+"&&e.dataset.parent==="+json.dumps(returned)+")")
    # Roving keyboard selection inspects without restoring document content.
    b.evaluate(node(branch)+'.focus()')
    key('ArrowDown')
    b.wait("document.activeElement?.dataset.eventId==="+json.dumps(returned))
    assert selected()==returned
    title('history_C')
    key('End')
    b.wait("document.activeElement?.dataset.eventId==="+json.dumps(root))
    assert selected()==root
    key('Home')
    b.wait("document.activeElement?.dataset.eventId==="+json.dumps(branch))
    inspect(returned)
    b.click('Restore selected event')
    title('history_A');assert current()==returned
    b.evaluate("document.querySelector('button[title=Undo]').click()")
    title('history_C');assert current()==branch
    b.evaluate("document.querySelector('button[title=Redo]').click()")
    title('history_A');assert current()==returned
    inspect(returned)
    b.click('Open selected event as new document')
    title('history_A')
    b.wait("document.querySelectorAll('.document-entry').length===3")
    fork=current();assert fork not in {root,middle,returned,branch}
    assert b.evaluate("document.querySelectorAll('.history-graph-event').length===1")
    assert b.evaluate("document.querySelector('button[title=Undo]').disabled")
    assert b.evaluate("Array.from(document.querySelectorAll('.history-graph-edge')).some(e=>e.dataset.child==="+json.dumps(fork)+"&&e.dataset.parent==="+json.dumps(returned)+"&&e.dataset.continuation==='document')")
    assert b.evaluate("document.querySelector('.history-graph-parent-list').textContent.includes('Parent in document [document-') && document.querySelector('.history-graph-parent-list').textContent.includes('history_A')")
    assert b.evaluate("document.querySelector('.history-graph-parent-list').textContent.includes("+json.dumps(returned)+")")
    # Exact event links persist across the real native archive autosave/reload path.
    b.wait("document.body.innerText.includes('Saved in this browser')")
    b.call('Page.reload')
    title('history_A')
    b.wait("document.body.innerText.includes('Saved in this browser')")
    assert current()==fork
    assert b.evaluate("document.querySelector('.history-graph-parent-list').textContent.includes("+json.dumps(returned)+")")
    # Reopen original document with the same displayed name; identity stays separate.
    b.evaluate("Array.from(document.querySelectorAll('.document-entry')).find(e=>e.querySelector('strong')?.textContent==='history_A').click()")
    assert current()==returned
    assert b.evaluate("document.querySelectorAll('.history-graph-event').length===4")
    # An exact parent event can be found without restoring the molecule.
    fill('input[aria-label="History event ID"]',root)
    b.click('Find event')
    assert selected()==root
    assert current()==returned
    # Exercise real page boundaries, not merely a CSS-limited first100-node prefix.
    for i in range(102):
        rename('history_step_'+str(i))
    last=current()
    assert b.evaluate("document.querySelectorAll('.history-graph-event').length===100")
    assert b.evaluate("Array.from(document.querySelectorAll('.history-graph-edge')).some(e=>e.dataset.continuation==='page')")
    b.click('Next 100 events')
    b.wait("document.querySelectorAll('.history-graph-event').length===6")
    assert b.evaluate("document.querySelectorAll('.history-graph-event[tabindex=\"0\"]').length===1")
    assert b.evaluate(node(root)+'!==null')
    assert not b.evaluate("document.querySelector('.history-panel').textContent.includes('on : click')")
    b.click('Show current event')
    b.wait("document.activeElement?.dataset.eventId==="+json.dumps(last))
    key('End')
    b.wait("document.activeElement?.dataset.eventId==="+json.dumps(root))
    assert selected()==root
    assert b.evaluate("document.querySelectorAll('.history-graph-event').length===6")
    # History graph scrolls inside its inspector at mobile width.
    if mobile_enabled():
        b.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True})
    b.wait("document.querySelector('.history-graph-scroll')")
    if mobile_enabled():
        assert b.evaluate('document.documentElement.scrollWidth<=390')
        assert b.evaluate("document.querySelector('.history-graph-scroll').getBoundingClientRect().right<=390")
    b.screenshot('/tmp/collagene2-history-graph.png')
    assert not b.exceptions,b.exceptions
    assert not b.console_errors,b.console_errors
    print('PASS: literal branch edges, equal-content distinctevents, exact undo/redo, keyboard inspection, sourcebound fork/reload, graph pagination and mobile bounds'.replace("mobile", "mobile checked" if mobile_enabled() else "mobile deferred"))


if __name__=='__main__':main()
