#!/usr/bin/env python3
"""Stage a bounded compiled preview, keeping old versioned workers for open tabs.

This copies only validated static snapshots. It grants no application filesystem
capability, executes no asset and deliberately drops obsolete main bundles.
"""
import argparse
import base64
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit
from serve_preview import MAX_FILES, MAX_TOTAL_BYTES, snapshot

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
    worker_hashes = {active.group(1)}
    worker_hashes.update(path.split("/")[2] for path in assets if WORKER.fullmatch(path))
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


def stage(build, previous, output, verify_assets=False, previous_active_only=False):
    assets = dict(snapshot(build))
    for old_root in previous:
        old_assets = snapshot(old_root)
        old_active = None
        if previous_active_only:
            old_active = validate_assets(old_assets)["active_worker"]
        for path, asset in old_assets.items():
            if not WORKER.fullmatch(path):
                continue
            if path in assets and assets[path][0] != asset[0]:
                raise ValueError("A versioned worker path has conflicting bytes.")
            if old_active is not None and path.split("/")[2] != old_active:
                continue
            assets.setdefault(path, asset)
        del old_assets
    total = sum(len(asset[0]) for asset in assets.values())
    if len(assets) > MAX_FILES or total > MAX_TOTAL_BYTES:
        raise ValueError("Merged preview exceeds immutable server limits.")
    if verify_assets:
        validate_assets(assets)
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
    parser.add_argument("--previous", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--previous-active-only", action="store_true", help="Keep only each previous snapshot's active worker, rather than its full historical worker collection")
    parser.add_argument("--verify-assets", action="store_true", help="Validate index references, integrity digests and complete versioned workers before staging")
    parser.add_argument("--manifest", type=Path, help="Write static verification evidence to a new file (requires --verify-assets)")
    args = parser.parse_args()
    if args.manifest and not args.verify_assets:
        parser.error("--manifest requires --verify-assets")
    if args.manifest and args.manifest.exists():
        parser.error("Manifest destination already exists")
    files, total = stage(args.build, args.previous, args.output, args.verify_assets, args.previous_active_only)
    if args.manifest:
        evidence = validate_assets(snapshot(args.output))
        evidence["root"] = str(args.output.resolve())
        evidence["previous_active_only"] = args.previous_active_only
        with args.manifest.open("x") as target:
            json.dump(evidence, target, indent=2)
            target.write("\n")
    print(f"Staged {files} immutable assets, {total} bytes in {args.output}")


if __name__ == "__main__":
    main()
