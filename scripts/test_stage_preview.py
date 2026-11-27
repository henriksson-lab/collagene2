"""Offline asset-staging checks; these never start a server or make HTTP requests."""
import base64
import hashlib
from pathlib import Path
import tempfile
import unittest

from stage_preview import stage, validate_assets
from serve_preview import snapshot

WASM = b"\x00asm\x01\x00\x00\x00"


def build(root, version="a" * 64):
    root.mkdir()
    (root / "main.js").write_bytes(b"export {};")
    (root / "main.wasm").write_bytes(WASM)
    (root / "style.css").write_bytes(b"body { color: black; }")
    digest = base64.b64encode(hashlib.sha384(b"export {};").digest()).decode()
    (root / "index.html").write_text(
        f'<meta name="collagene-worker" content="/workers/{version}/analysis-worker_loader.js">'
        f'<script src="/main.js" integrity="sha384-{digest}"></script>'
        '<link rel="stylesheet" href="style.css?v=1">'
        '<link rel="preload" href="/main.wasm" as="fetch">'
    )
    worker = root / "workers" / version
    worker.mkdir(parents=True)
    (worker / "analysis-worker_loader.js").write_bytes(b"loader")
    (worker / "analysis-worker.js").write_bytes(b"worker")
    (worker / "analysis-worker_bg.wasm").write_bytes(WASM)


class PreviewStaging(unittest.TestCase):
    def test_current_assets_survive_staging_with_verified_references(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            fresh, output = parent / "fresh", parent / "output"
            build(fresh)
            stage(fresh, output)
            assets = snapshot(output)
            evidence = validate_assets(assets)
            self.assertEqual(evidence["active_worker"], "a" * 64)
            self.assertEqual(evidence["worker_hashes"], ["a" * 64])
            self.assertEqual(assets["/main.js"][0], b"export {};")
            self.assertFalse(evidence["application_browser_checks"])
            self.assertFalse(evidence["http_smoke_checks"])

    def test_invalid_bundle_is_rejected_before_output_is_created(self):
        for problem in ("integrity", "missing", "worker", "wasm", "namespace", "extra_worker"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory() as directory:
                parent = Path(directory)
                fresh, output = parent / "fresh", parent / "output"
                build(fresh)
                if problem == "integrity":
                    (fresh / "main.js").write_bytes(b"wrong bundle bytes")
                elif problem == "missing":
                    (fresh / "style.css").unlink()
                elif problem == "worker":
                    (fresh / "workers" / ("a" * 64) / "analysis-worker.js").unlink()
                elif problem == "wasm":
                    (fresh / "main.wasm").write_bytes(b"not wasm")
                elif problem == "extra_worker":
                    worker = fresh / "workers" / ("b" * 64)
                    worker.mkdir()
                    (worker / "analysis-worker_loader.js").write_bytes(b"obsolete loader")
                    (worker / "analysis-worker.js").write_bytes(b"obsolete worker")
                    (worker / "analysis-worker_bg.wasm").write_bytes(WASM)
                else:
                    p = fresh / "index.html"
                    p.write_text(p.read_text().replace("style.css?v=1", "%2e%2e/style.css"))
                with self.assertRaises(ValueError):
                    stage(fresh, output)
                self.assertFalse(output.exists())

    def test_existing_snapshot_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            fresh, output = parent / "fresh", parent / "output"
            build(fresh)
            output.mkdir()
            (output / "index.html").write_bytes(b"currently published")
            with self.assertRaises(FileExistsError):
                stage(fresh, output)
            self.assertEqual((output / "index.html").read_bytes(), b"currently published")


if __name__ == "__main__":
    unittest.main()
