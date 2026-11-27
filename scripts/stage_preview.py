#!/usr/bin/env python3
"""Stage one bounded current compiled preview.

This copies only validated static snapshots. It grants no application filesystem
capability, executes no asset.
"""
import argparse
import base64
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit
from serve_preview import snapshot

WORKER = re.compile(r"/workers/[0-9a-f]{64}/analysis-worker(?:_loader\.js|\.js|_bg\.wasm)\Z")


class AssetIndex(HTMLParser):
    def __init__(self):
        super().__init__()
        self.workers = []
        self.references = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "meta" and attrs.get("name") == "collagene-worker":
            self.workers.append(attrs.get("content", ""))
        if tag == "script" and attrs.get("src"):
            self.references.append((attrs["src"], attrs.get("integrity")))
        if (tag == "link" and attrs.get("href")
                and set(attrs.get("rel", "").split()) & {"stylesheet", "modulepreload", "preload"}):
            self.references.append((attrs["href"], attrs.get("integrity")))


def asset_path(reference):
    url = urlsplit(reference)
    if url.scheme or url.netloc:
        raise ValueError("Compiled preview references must use local snapshot assets.")
    path = "/" + unquote(url.path, errors="strict").lstrip("/")
    if any(part in (".", "..") or part.startswith(".") for part in path.split("/") if part):
        raise ValueError("Compiled preview reference is outside the served asset namespace.")
    return path


def validate_assets(assets):
    """Check static references and declared hashes; never execute or request assets."""
    index = AssetIndex()
    index.feed(assets["/index.html"][0].decode("utf-8"))
    if len(index.workers) != 1:
        raise ValueError("Expected one active analysis-worker meta reference.")
    active = re.fullmatch(r"/workers/([0-9a-f]{64})/analysis-worker_loader\.js", index.workers[0])
    if active is None:
        raise ValueError("Active worker must have a versioned snapshot path.")
    if any(path.startswith("/workers/") and not WORKER.fullmatch(path) for path in assets):
        raise ValueError("Unexpected asset in the analysis-worker namespace.")
    worker_hashes = {path.split("/")[2] for path in assets if WORKER.fullmatch(path)}
    if worker_hashes != {active.group(1)}:
        raise ValueError("A preview must contain only its current analysis worker.")
    for version in worker_hashes:
        for name in ("analysis-worker_loader.js", "analysis-worker.js", "analysis-worker_bg.wasm"):
            if f"/workers/{version}/{name}" not in assets:
                raise ValueError("Versioned analysis-worker files are incomplete.")
    for reference, integrity in index.references:
        path = asset_path(reference)
        if path not in assets:
            raise ValueError(f"Referenced preview asset is missing: {path}")
        for token in (integrity or "").split():
            algorithm, separator, encoded = token.partition("-")
            if not separator or algorithm not in ("sha256", "sha384", "sha512"):
                raise ValueError(f"Unsupported asset integrity digest: {path}")
            try:
                declared = base64.b64decode(encoded, validate=True)
            except ValueError as error:
                raise ValueError(f"Malformed asset integrity digest: {path}") from error
            if hashlib.new(algorithm, assets[path][0]).digest() != declared:
                raise ValueError(f"Asset integrity mismatch: {path}")
    for path, (data, _, _) in assets.items():
        if path.endswith(".wasm") and data[:8] != b"\x00asm\x01\x00\x00\x00":
            raise ValueError(f"Invalid WebAssembly header: {path}")
    return {
        "assets": len(assets),
        "bytes": sum(len(asset[0]) for asset in assets.values()),
        "index_sha256": hashlib.sha256(assets["/index.html"][0]).hexdigest(),
        "active_worker": active.group(1),
        "worker_hashes": sorted(worker_hashes),
        "static_integrity_checks": True,
        "application_browser_checks": False,
        "http_smoke_checks": False,
    }


def stage(build, output):
    assets = snapshot(build)
    validate_assets(assets)
    total = sum(len(asset[0]) for asset in assets.values())
    # A new destination prevents accidental overwrites of an existing release.
    output = Path(output)
    output.mkdir(parents=False, exist_ok=False)
    for path, (data, _, _) in sorted(assets.items()):
        destination = output / path.lstrip("/")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as target:
            target.write(data)
    return len(assets), total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, help="Write static verification evidence to a new file")
    args = parser.parse_args()
    if args.manifest and args.manifest.exists():
        parser.error("Manifest destination already exists")
    files, total = stage(args.build, args.output)
    if args.manifest:
        evidence = validate_assets(snapshot(args.output))
        evidence["root"] = str(args.output.resolve())
        with args.manifest.open("x") as target:
            json.dump(evidence, target, indent=2)
            target.write("\n")
    print(f"Staged {files} immutable assets, {total} bytes in {args.output}")


if __name__ == "__main__":
    main()
