#!/usr/bin/env python3
"""Check embedded stdio mode of a built desktop binary; never start its GUI."""
import argparse
import json
import queue
import subprocess
import threading
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    args = parser.parse_args()
    process = subprocess.Popen(
        [str(args.binary.resolve()), "--mcp"],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, encoding="utf-8", bufsize=1,
    )
    replies = queue.Queue()
    diagnostics = bytearray()

    def read_stdout():
        try:
            limit = 64 * 1024 * 1024
            while line := process.stdout.readline(limit + 1):
                if len(line) > limit:
                    raise ValueError("embedded MCP output exceeded the check limit")
                replies.put(json.loads(line))
        except Exception as error:
            replies.put(error)
        finally:
            replies.put(None)

    def read_stderr():
        while chunk := process.stderr.read(4096):
            raw = chunk.encode("utf-8", errors="replace")
            diagnostics.extend(raw[: max(0, 32768 - len(diagnostics))])

    threading.Thread(target=read_stdout, daemon=True).start()
    threading.Thread(target=read_stderr, daemon=True).start()
    next_id = 0

    def send(method, params, notification=False):
        nonlocal next_id
        next_id += 1
        request = {"jsonrpc": "2.0", "method": method, "params": params}
        if not notification:
            request["id"] = next_id
        process.stdin.write(json.dumps(request) + "\n")
        process.stdin.flush()
        if notification:
            return None
        while True:
            response = replies.get(timeout=25)
            if response is None or isinstance(response, Exception):
                raise RuntimeError("desktop --mcp did not produce valid stdio MCP")
            if response.get("id") == next_id:
                if "error" in response:
                    raise RuntimeError(f"embedded MCP rejected {method}: {response['error']}")
                return response["result"]

    def tool(name, arguments):
        result = send("tools/call", {"name": name, "arguments": arguments})
        if result.get("isError"):
            raise RuntimeError(f"embedded MCP tool failed: {name}")
        return result["structuredContent"]

    try:
        initialized = send("initialize", {
            "protocolVersion": "2025-11-25", "capabilities": {},
            "clientInfo": {"name": "desktop-embedded-check", "version": "1"},
        })
        assert initialized["serverInfo"]["name"] == "collagene"
        send("notifications/initialized", {}, notification=True)
        available = {value["name"] for value in send("tools/list", {})["tools"]}
        assert {"import_document", "verify_plan", "inspect_restriction"} <= available
        document = tool("import_document", {"text": ">embedded-source\nAAAAGAATTCAAAA\n"})["documents"][0]
        digest = tool("inspect_restriction", {
            "document": {"id": document["id"], "revision": document["revision"]},
            "source_role": "embedded-source", "enzymes": ["EcoRI"],
        })
        assert digest["scan_complete"] and digest["products_determined"]
        assert digest["site_count"] == 1 and digest["fragment_count"] == 2
        process.stdin.close()
        assert process.wait(timeout=15) == 0
        print("Embedded desktop MCP: stdio initialization, cloning tools and virtual digest passed; GUI not launched.")
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=10)
        if process.returncode:
            print(diagnostics.decode("utf-8", errors="replace"))


if __name__ == "__main__":
    main()
