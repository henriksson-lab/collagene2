#!/usr/bin/env python3
"""Read-only static GUI preview from an explicit immutable asset snapshot.

No application filesystem API, file uploads, directory listing or disk writes.
The browser's explicitly selected documents remain in that visitor's browser.
"""
import argparse
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import mimetypes
import os
from pathlib import Path
import stat
from types import MappingProxyType
from urllib.parse import unquote, urlsplit

MAX_FILES = 4096
MAX_FILE_BYTES = 128 * 1024 * 1024
# Bounded immutable assets, including prior versioned workers for open tabs.
MAX_TOTAL_BYTES = 385 * 1024 * 1024


def snapshot(directory):
    root = Path(directory).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("The asset root must be an existing directory.")
    assets = {}
    total = 0
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if any(part.startswith(".") for part in relative.parts):
            continue
        info = path.lstat()
        if stat.S_ISLNK(info.st_mode):
            raise ValueError("Asset symlinks are not permitted.")
        if stat.S_ISDIR(info.st_mode):
            continue
        if not stat.S_ISREG(info.st_mode) or path.resolve(strict=True) != path:
            raise ValueError("Only regular files inside the asset root are permitted.")
        if len(assets) >= MAX_FILES or info.st_size > MAX_FILE_BYTES:
            raise ValueError("Preview asset limit exceeded.")
        flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_BINARY", 0)
        with os.fdopen(os.open(path, flags), "rb") as source:
            opened = os.fstat(source.fileno())
            if not stat.S_ISREG(opened.st_mode) or (opened.st_dev, opened.st_ino) != (info.st_dev, info.st_ino):
                raise ValueError("Asset changed while opening the snapshot.")
            budget = min(MAX_FILE_BYTES, MAX_TOTAL_BYTES - total)
            data = source.read(budget + 1)
            if len(data) > budget:
                raise ValueError("Preview asset limit exceeded.")
        total += len(data)
        mime = "application/wasm" if path.suffix == ".wasm" else mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        assets["/" + relative.as_posix()] = (data, mime, '"' + hashlib.sha256(data).hexdigest() + '"')
    if "/index.html" not in assets:
        raise ValueError("The asset root must contain index.html.")
    return MappingProxyType(assets)


def handler_for(assets):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *_):
            # Do not retain visitor URLs/query strings or document names.
            pass

        def serve(self, body):
            try:
                path = unquote(urlsplit(self.path).path, errors="strict")
            except (ValueError, UnicodeError):
                self.send_error(400)
                return
            if any(part in (".", "..") or part.startswith(".") for part in path.split("/") if part):
                self.send_error(404)
                return
            asset = assets.get("/index.html" if path == "/" else path)
            if asset is None:
                self.send_error(404)
                return
            data, mime, etag = asset
            self.send_response(200)
            self.send_header("Content-Type", mime)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("ETag", etag)
            self.send_header("Cache-Control", "no-cache")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            if body:
                self.wfile.write(data)

        def do_GET(self):
            self.serve(True)

        def do_HEAD(self):
            self.serve(False)

        def do_POST(self):
            self.send_error(405, "This preview is read-only")

        do_PUT = do_POST
        do_PATCH = do_POST
        do_DELETE = do_POST

    return Handler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Only this compiled asset directory is read at startup; writes are disabled")
    parser.add_argument("--bind", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8082)
    args = parser.parse_args()
    assets = snapshot(args.root)
    server = ThreadingHTTPServer((args.bind, args.port), handler_for(assets))
    print(f"Read-only GUI preview on {args.bind}:{args.port}; {len(assets)} snapshotted assets; no filesystem API", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
