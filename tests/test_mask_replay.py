import sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import verify_recovered_mask as mask

class MaskReplayChecks(unittest.TestCase):
    def test_merge_preserves_union_not_sum(self):
        self.assertEqual(mask.merge([(0,5),(3,8),(8,10),(12,14)]), [(0,10),(12,14)])
    def test_intersection_uses_half_open_intervals(self):
        self.assertEqual(mask.intersect({'x':[(0,5),(8,12)]},{'x':[(5,10)]}), {'x':[(8,10)]})
    def test_missing_input_fails_without_output(self):
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/'output'
            with self.assertRaises(FileNotFoundError):mask.verify_recovered_mask(Path(t)/'missing',out)
            self.assertFalse(out.exists())
    def test_wrong_contract_hash_fails_without_output(self):
        with tempfile.TemporaryDirectory() as t:
            data=Path(t)/'data';(data/'provenance').mkdir(parents=True)
            (data/'provenance/SD07_contract.json').write_text('{}')
            out=Path(t)/'output'
            with self.assertRaisesRegex(ValueError,'identity mismatch'):mask.verify_recovered_mask(data,out)
            self.assertFalse(out.exists())
if __name__=='__main__':unittest.main()
