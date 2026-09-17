import unittest, tempfile, json
from pathlib import Path
import numpy as np
from .predictions import score_predictions

class FrozenPredictionChecks(unittest.TestCase):
    def fixture(self,root,pred=None,observed=None,pt=None,ot=None,pf=None,of=None):
        pred=np.eye(3) if pred is None else pred;observed=np.eye(3) if observed is None else observed
        np.save(root/'pred.npy',pred);np.save(root/'obs.npy',observed)
        cfg={'dataset_id':'test','gallery_targets':['a','b','c'],'evaluation_provenance':{'scope':'unit check'},'predictions':{'values':'pred.npy','targets':pt or ['a','b','c'],'features':pf or ['x','y','z'],'output_space':'same','readout':'delta','provenance':{'scope':'unit check'}},'observed':{'values':'obs.npy','targets':ot or ['a','b','c'],'features':of or ['x','y','z'],'output_space':'same','readout':'delta','provenance':{'scope':'unit check'}}}
        cp=root/'config.json';cp.write_text(json.dumps(cfg));return cp
    def test_reordered_axes_recover_identity(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);cp=self.fixture(r,observed=np.eye(3)[[2,0,1]][:,[1,2,0]],ot=['c','a','b'],of=['y','z','x'])
            report=score_predictions(cp,r/'out');self.assertEqual(report['summary']['mean_correct_target_rank'],1)
            np.testing.assert_array_equal(np.load(r/'out/score_matrix.npz')['scores'],np.eye(3))
    def test_missing_prediction_retains_gallery_and_coverage(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);cp=self.fixture(r,pred=np.eye(3)[[0,2]],pt=['a','c']);report=score_predictions(cp,r/'out')
            self.assertEqual(report['fixed_gallery_size'],3);self.assertEqual(report['missing_prediction_targets'],['b']);self.assertEqual(report['summary']['n_finite_query_rows'],2)
            self.assertTrue(np.isnan(np.load(r/'out/score_matrix.npz')['scores'][1]).all())
    def test_coordinate_mismatch_and_duplicate_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);cp=self.fixture(r,of=['x','y','unmatched'])
            with self.assertRaisesRegex(ValueError,'intersection'):score_predictions(cp,r/'out')
            self.assertFalse((r/'out').exists());cp=self.fixture(r,pt=['a','a','c'])
            with self.assertRaisesRegex(ValueError,'duplicate'):score_predictions(cp,r/'out')
    def test_zero_direction_not_fabricated_as_chance(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);cp=self.fixture(r,pred=np.zeros((3,3)));report=score_predictions(cp,r/'out')
            self.assertEqual(report['summary']['n_finite_query_rows'],0);self.assertIsNone(report['summary']['mean_correct_target_rank'])
    def test_no_predictions_is_explicit_unavailable_not_dropped(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);cp=self.fixture(r,pred=np.zeros((0,3)));c=json.loads(cp.read_text());c['predictions']['targets']=[];cp.write_text(json.dumps(c))
            report=score_predictions(cp,r/'out');self.assertEqual(report['missing_prediction_targets'],['a','b','c']);self.assertEqual(report['summary']['n_finite_query_rows'],0)
    def test_nonfinite_observed_invalidates_fixed_gallery(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);obs=np.eye(3);obs[1,0]=np.nan;cp=self.fixture(r,observed=obs);report=score_predictions(cp,r/'out')
            self.assertEqual(report['summary']['n_finite_query_rows'],0)
    def test_output_space_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);cp=self.fixture(r);c=json.loads(cp.read_text());c['observed']['output_space']='different';cp.write_text(json.dumps(c))
            with self.assertRaisesRegex(ValueError,'output_space'):score_predictions(cp,r/'out')

if __name__=='__main__':unittest.main()
