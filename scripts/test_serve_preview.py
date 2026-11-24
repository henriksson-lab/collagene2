"""Development checks for the public read-only preview boundary."""
from contextlib import contextmanager
from http.client import HTTPConnection
from http.server import ThreadingHTTPServer
from pathlib import Path
import tempfile
import threading
import unittest

from serve_preview import handler_for, snapshot


@contextmanager
def running(assets):
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler_for(assets))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server.server_port
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


class PreviewBoundary(unittest.TestCase):
    def test_requests_cannot_escape_list_or_write(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            root = parent / "public"
            root.mkdir()
            (root / "index.html").write_text("GUI")
            (parent / "outside.txt").write_text("outside sentinel")
            (root / ".private").write_text("private sentinel")
            assets = snapshot(root)
            with running(assets) as port:
                for method, path, expected in [
                    ("GET", "/", 200),
                    ("HEAD", "/index.html", 200),
                    ("GET", "/../outside.txt", 404),
                    ("GET", "/%2e%2e/outside.txt", 404),
                    ("GET", "/.private", 404),
                    ("GET", "/missing/", 404),
                    ("POST", "/index.html", 405),
                    ("PUT", "/new.txt", 405),
                    ("DELETE", "/index.html", 405),
                ]:
                    with self.subTest(method=method, path=path):
                        connection = HTTPConnection("127.0.0.1", port)
                        connection.request(method, path)
                        response = connection.getresponse()
                        body = response.read()
                        self.assertEqual(response.status, expected)
                        self.assertNotIn(b"outside sentinel", body)
                        self.assertNotIn(b"private sentinel", body)
                        connection.close()
            self.assertEqual((root / "index.html").read_text(), "GUI")
            self.assertFalse((root / "new.txt").exists())

    def test_snapshot_stays_read_only_when_disk_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "index.html"
            path.write_text("original")
            assets = snapshot(root)
            path.write_text("changed")
            self.assertEqual(assets["/index.html"][0], b"original")
            with self.assertRaises(TypeError):
                assets["/new"] = b"injected"

    def test_symlinks_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            root = parent / "public"
            root.mkdir()
            (root / "index.html").write_text("GUI")
            outside = parent / "outside.txt"
            outside.write_text("sentinel")
            try:
                (root / "leak.txt").symlink_to(outside)
            except OSError:
                self.skipTest("OS account cannot create symlinks")
            with self.assertRaises(ValueError):
                snapshot(root)


if __name__ == "__main__":
    unittest.main()
