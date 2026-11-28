#!/usr/bin/env python3
"""Check actual headless desktop payloads and stage the requested CI artifacts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tarfile


MAX_FILE_BYTES = 1024 * 1024 * 1024
MAX_PLATFORM_BYTES = 2 * MAX_FILE_BYTES


def one(pattern):
    paths = list(Path("target").glob(pattern))
    if len(paths) != 1:
        raise SystemExit("missing or ambiguous desktop payload: " + pattern)
    return paths[0]


def windows_dependencies(binary):
    installer = Path(os.environ["ProgramFiles(x86)"]) / "Microsoft Visual Studio/Installer/vswhere.exe"
    installation = subprocess.check_output(
        [str(installer), "-latest", "-products", "*", "-requires",
         "Microsoft.VisualStudio.Component.VC.Tools.x86.x64", "-property", "installationPath"],
        text=True,
    ).strip()
    if not installation:
        raise SystemExit("Visual C++ toolchain installation was not found")
    tools = Path(installation) / "VC/Tools/MSVC"
    candidates = [p for p in tools.iterdir() if re.fullmatch(r"\d+(?:\.\d+)+", p.name)]
    if not candidates:
        raise SystemExit("Visual C++ toolchain version was not found")
    version = max(candidates, key=lambda p: tuple(map(int, p.name.split("."))))
    dumpbin = version / "bin/Hostx64/x64/dumpbin.exe"
    report = subprocess.check_output(
        [str(dumpbin), "/DEPENDENTS", str(binary)], text=True, errors="replace",
    )
    dependencies = sorted(set(re.findall(r"\b[\w.-]+\.dll\b", report, re.IGNORECASE)))
    if not dependencies:
        raise SystemExit("could not inspect portable Windows executable dependencies")
    system_dlls = Path(os.environ["SystemRoot"]) / "System32"
    unavailable = [name for name in dependencies if re.fullmatch(
        r"(?:vcruntime\d.*|msvcp\d.*|concrt\d.*|ucrtbased|msvcr\d.*|webview2loader)\.dll",
        name, re.IGNORECASE,
    ) or not (name.lower().startswith(("api-ms-win-", "ext-ms-win-"))
              or (system_dlls / name).is_file())]
    if unavailable:
        raise SystemExit("portable executable has non-system runtime imports: " + ", ".join(unavailable))
    return dependencies


def copy_payload(source, output):
    if source.is_symlink() or not source.is_file() or source.stat().st_size == 0:
        raise SystemExit("package payload must be a nonempty regular file: " + str(source))
    if source.stat().st_size > MAX_FILE_BYTES:
        raise SystemExit("package exceeds 1 GiB artifact limit")
    target = output / source.name
    if target.exists():
        raise SystemExit("duplicate output package: " + source.name)
    shutil.copy2(source, target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=["development", "release"], required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--frontend-artifact", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.source_revision):
        raise SystemExit("source revision must be the exact Git commit")
    if not args.frontend_artifact or len(args.frontend_artifact) > 256:
        raise SystemExit("frontend artifact identity is required")
    profile = "debug" if args.profile == "development" else "release"
    system = platform.system()
    output = Path("desktop-packages")
    # A restored native cache can contain old bundle outputs. The workflow deletes
    # its selected bundle directory before building; staging itself never merges.
    output.mkdir()
    prerequisites = []
    dependencies = []
    if system == "Windows":
        target = "x86_64-pc-windows-msvc"
        binary = one(f"{target}/{profile}/collagene-desktop.exe")
        dependencies = windows_dependencies(binary)
        prerequisites = ["Microsoft Edge WebView2 Runtime (installed in Windows)"]
        payloads = [binary]
    elif system == "Darwin":
        apps = list((Path("target") / profile / "bundle/macos").glob("*.app"))
        if len(apps) != 1 or not apps[0].is_dir() or apps[0].is_symlink():
            raise SystemExit("missing or ambiguous macOS application bundle")
        app = apps[0]
        binary = app / "Contents/MacOS/collagene-desktop"
        payloads = [one(f"{profile}/bundle/dmg/*.dmg")]
        with tarfile.open(output / (app.name + ".tar.gz"), "w:gz") as archive:
            # Preserve application executable bits and internal symlinks.
            archive.add(app, arcname=app.name)
    elif system == "Linux":
        package = one(f"{profile}/bundle/deb/*.deb")
        extracted = Path(os.environ["RUNNER_TEMP"]) / "collagene-deb-payload"
        subprocess.run(["dpkg-deb", "--extract", str(package), str(extracted)], check=True)
        binary = extracted / "usr/bin/collagene-desktop"
        payloads = [package]
        prerequisites = ["Runtime packages declared by the Debian package"]
    else:
        raise SystemExit("unsupported desktop packaging platform")
    if binary.is_symlink() or not binary.is_file():
        raise SystemExit("missing actual desktop executable")
    subprocess.run([sys.executable, "scripts/desktop_mcp_smoke.py", str(binary)], check=True)
    for payload in payloads:
        copy_payload(payload, output)
    expected = 2 if system == "Darwin" else 1
    files = list(output.iterdir())
    if len(files) != expected:
        raise SystemExit("unexpected staged platform artifact count")
    checksums = {}
    total = 0
    for path in files:
        size = path.stat().st_size
        if size == 0 or size > MAX_FILE_BYTES:
            raise SystemExit("package is empty or exceeds 1 GiB")
        total += size
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        checksums[path.name] = digest.hexdigest()
    if total > MAX_PLATFORM_BYTES:
        raise SystemExit("package artifacts exceed 2 GiB total")
    metadata = {
        "source_revision": args.source_revision,
        "frontend_artifact": args.frontend_artifact,
        "profile": args.profile,
        "os": system,
        "architecture": platform.machine(),
        "cargo_tauri": "2.11.3",
        "signed": False,
        "runtime_prerequisites": prerequisites,
        "windows_imported_system_dlls": dependencies,
        "sha256": checksums,
    }
    (output / "BUILD.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Staged {expected} requested desktop artifacts ({total} bytes) for {args.source_revision}")


if __name__ == "__main__":
    main()
