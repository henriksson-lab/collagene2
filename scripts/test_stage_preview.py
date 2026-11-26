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
    def test_fresh_assets_and_exact_old_workers_survive_staging(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            fresh, old, output = (parent / name for name in ("fresh", "old", "output"))
            build(fresh)
            build(old, "b" * 64)
            (old / "main.js").write_bytes(b"obsolete main")
            stage(fresh, [old], output, verify_assets=True)
            assets = snapshot(output)
            evidence = validate_assets(assets)
            self.assertEqual(evidence["active_worker"], "a" * 64)
            self.assertEqual(evidence["worker_hashes"], ["a" * 64, "b" * 64])
            self.assertEqual(assets["/main.js"][0], b"export {};")
            for path, asset in snapshot(old).items():
                if path.startswith("/workers/"):
                    self.assertEqual(assets[path][0], asset[0])
            self.assertFalse(evidence["application_browser_checks"])
            self.assertFalse(evidence["http_smoke_checks"])

    def test_bounded_retention_keeps_each_previous_active_worker_without_copying_ancient_versions(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            fresh, old, ancient, old_merged, output = (parent / name for name in ("fresh", "old", "ancient", "old_merged", "output"))
            build(fresh)
            build(old, "b" * 64)
            build(ancient, "c" * 64)
            stage(old, [ancient], old_merged, verify_assets=True)
            stage(fresh, [old_merged], output, verify_assets=True, previous_active_only=True)
            assets = snapshot(output)
            self.assertEqual(validate_assets(assets)["worker_hashes"], ["a" * 64, "b" * 64])
            self.assertEqual(assets["/workers/" + "b" * 64 + "/analysis-worker.js"][0], b"worker")
            self.assertFalse(any(path.startswith("/workers/" + "c" * 64 + "/") for path in assets))

    def test_invalid_bundle_is_rejected_before_output_is_created(self):
        for problem in ("integrity", "missing", "worker", "wasm", "namespace"):
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
                else:
                    p = fresh / "index.html"
                    p.write_text(p.read_text().replace("style.css?v=1", "%2e%2e/style.css"))
                with self.assertRaises(ValueError):
                    stage(fresh, [], output, verify_assets=True)
                self.assertFalse(output.exists())

    def test_same_worker_identity_cannot_replace_old_tab_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            fresh, old, output = (parent / name for name in ("fresh", "old", "output"))
            build(fresh)
            build(old)
            (old / "workers" / ("a" * 64) / "analysis-worker.js").write_bytes(b"conflicting worker")
            with self.assertRaises(ValueError):
                stage(fresh, [old], output, verify_assets=True)
            self.assertFalse(output.exists())

    def test_existing_snapshot_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            fresh, output = parent / "fresh", parent / "output"
            build(fresh)
            output.mkdir()
            (output / "index.html").write_bytes(b"currently published")
            with self.assertRaises(FileExistsError):
                stage(fresh, [], output, verify_assets=True)
            self.assertEqual((output / "index.html").read_bytes(), b"currently published")


if __name__ == "__main__":
    unittest.main()
