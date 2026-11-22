import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import summarize_resource as s
import verify_sd11 as v
class Semantics(unittest.TestCase):
    def test_partition_boundaries(self):
        self.assertEqual(s.classify(127,1,True),'SELECTED')
        self.assertEqual(s.classify(126,1,True),'LOW_CALL_COUNT_ONLY')
        self.assertEqual(s.classify(127,1,False),'OUTSIDE_MASK_ONLY')
        self.assertEqual(s.classify(126,1,False),'BOTH')
    def test_no_observations_is_not_reference(self):
        self.assertEqual(s.classify(0,0,True),'NO_RETAINED')
        self.assertEqual(s.classify(1,0,False),'MONOMORPHIC_REF')
    def test_variable_is_retained(self):
        data=[dict(fsp='example',host=h,reaction='V') for h in s.HOSTS]
        for r in s.reaction_counts(data,{'example'}):self.assertEqual((r['R'],r['S'],r['V']),(0,0,1))
    def test_low_depth_is_unresolved(self):
        d={k:0 for k in v.COUNT+v.OTHER};d.update(depth=4,ALT_F=4)
        self.assertEqual(v.label(d,'1')['supported_allele'],'.')
    def test_deletions_are_in_denominator(self):
        d={k:0 for k in v.COUNT+v.OTHER};d.update(depth=10,ALT_F=8,DEL=2)
        self.assertEqual(v.label(d,'1')['support_status'],'MIXED_OR_NONBASE')
    def test_missing_stays_missing_even_when_supported(self):
        d={k:0 for k in v.COUNT+v.OTHER};d.update(depth=10,REF_F=9,ALT_R=1)
        r=v.label(d,'.');self.assertEqual(r['baseline_relation'],'BASELINE_MISSING')
        self.assertEqual(r['supported_allele'],'0')
    def test_pileup_nonbase_tokens(self):
        d=v.unpack('scaffold\t1\tA\t5\t.,Gg*\tIIIII\tIIIII\n',{'ref':'A','alt':'G'})
        self.assertEqual([d[k] for k in ['REF_F','REF_R','ALT_F','ALT_R','DEL']],[1]*5)
    def test_pileup_truncation_fails(self):
        with self.assertRaises(AssertionError):v.unpack('s\t1\tA\t5\t.,Gg\tIIIII\tIIIII\n',{'ref':'A','alt':'G'})
if __name__=='__main__':unittest.main()
