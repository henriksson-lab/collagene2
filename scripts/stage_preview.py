#!/usr/bin/env python3
"""Stage a bounded compiled preview, keeping old versioned workers for open tabs.

This copies only validated static snapshots. It grants no application filesystem
capability, executes no asset and deliberately drops obsolete main bundles.
"""
import argparse
from pathlib import Path
import re
from serve_preview import MAX_FILES, MAX_TOTAL_BYTES, snapshot

WORKER = re.compile(r"/workers/[0-9a-f]{64}/analysis-worker(?:_loader\.js|\.js|_bg\.wasm)\Z")


def stage(build, previous, output):
    assets = dict(snapshot(build))
    for old_root in previous:
        old_assets = snapshot(old_root)
        for path, asset in old_assets.items():
            if not WORKER.fullmatch(path):
                continue
            if path in assets and assets[path][0] != asset[0]:
                raise ValueError("A versioned worker path has conflicting bytes.")
            assets.setdefault(path, asset)
        del old_assets
    total = sum(len(asset[0]) for asset in assets.values())
    if len(assets) > MAX_FILES or total > MAX_TOTAL_BYTES:
        raise ValueError("Merged preview exceeds immutable server limits.")
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
    args = parser.parse_args()
    files, total = stage(args.build, args.previous, args.output)
    print(f"Staged {files} immutable assets, {total} bytes in {args.output}")


if __name__ == "__main__":
    main()
