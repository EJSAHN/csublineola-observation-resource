import hashlib,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import replay_sd13 as r
class TargetedReplayChecks(unittest.TestCase):
    def fixture(self,root):
        (root/'data.txt').write_text('frozen')
        (root/'SHA256SUMS.txt').write_text(hashlib.sha256(b'frozen').hexdigest()+'  data.txt\n')
    def test_missing_input_fails(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaises(FileNotFoundError):r.replay(Path(t)/'missing',Path(t)/'out')
    def test_changed_file_fails(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.fixture(p);(p/'data.txt').write_text('changed')
            with self.assertRaises(ValueError):r.check_manifest(p,'SHA256SUMS.txt')
    def test_generated_cache_ignored(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.fixture(p);(p/'__pycache__').mkdir();(p/'__pycache__/example.pyc').write_bytes(b'local')
            self.assertEqual(r.check_manifest(p,'SHA256SUMS.txt'),1)
    def test_unlisted_source_fails(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.fixture(p);(p/'unexpected.py').write_text('pass')
            with self.assertRaises(ValueError):r.check_manifest(p,'SHA256SUMS.txt')
    def test_path_escape_fails(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);(p/'SHA256SUMS.txt').write_text('0'*64+'  ../escaped\n')
            with self.assertRaises(ValueError):r.check_manifest(p,'SHA256SUMS.txt')
    def test_cache_in_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);(p/'__pycache__').mkdir();(p/'__pycache__/example.pyc').write_bytes(b'local');(p/'SHA256SUMS.txt').write_text(hashlib.sha256(b'local').hexdigest()+'  __pycache__/example.pyc\n')
            with self.assertRaises(ValueError):r.check_manifest(p,'SHA256SUMS.txt')
if __name__=='__main__':unittest.main()
