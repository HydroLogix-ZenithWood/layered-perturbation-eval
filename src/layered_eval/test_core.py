import unittest
import numpy as np
from .core import *

class ScientificKernelChecks(unittest.TestCase):
    def test_identity_reordering(self):
        x=np.array([[20],[10],[30]])
        np.testing.assert_array_equal(align_rows(x,['b','a','c'],['a','b','c']).ravel(),[10,20,30])
        with self.assertRaises(ValueError):align_rows(x,['a','a','c'],['a','c'])
        with self.assertRaises(ValueError):align_rows(x,['b','a','c'],['d'])
    def test_ties_and_correct_map(self):
        m=summarize_scores(np.ones((5,5)));self.assertEqual(m['mean_correct_target_rank'],.5);self.assertEqual(m['top1_fraction_tie_adjusted'],.2)
        self.assertEqual(summarize_scores(np.eye(5))['mean_correct_target_rank'],1)
        self.assertEqual(summarize_scores(np.eye(5)[:,::-1],np.arange(4,-1,-1))['mean_correct_target_rank'],1)
    def test_zero_nan(self):
        s=cosine_matrix(np.zeros((3,2)),np.ones((3,2)))
        self.assertTrue(np.isnan(s).all());self.assertFalse(summarize_scores(s)['complete']);self.assertEqual(summarize_scores(s)['n_finite_query_rows'],0)
        with self.assertRaises(ValueError):conditional_label_null(s)
        y=np.eye(3);y[1,1]=np.nan;self.assertTrue(np.isnan(cosine_matrix(y,np.eye(3))[1]).all())
    def test_leakage(self):
        assert_disjoint([1,2],[3,4],[5])
        with self.assertRaises(ValueError):assert_disjoint([1,2],[2,3])
        with self.assertRaises(ValueError):assert_disjoint([1,1])
    def test_equal_guide_weights_and_gram(self):
        rng=np.random.default_rng(4);x=rng.normal(size=(40,13))
        a={'g1':{'a':[0,1],'b':[2,3,4]},'g2':{'c':[5,6]}}
        b={'g1':{'d':[7,8]},'g2':{'e':[9,10],'f':[11,12,13]}}
        ca={'c1':[20,21],'c2':[22,23,24]};cb={'c1':[25,26],'c2':[27,28,29]}
        u=aggregate_responses(x,a,ca);v=aggregate_responses(x,b,cb)
        wa=response_weights(40,a,ca);wb=response_weights(40,b,cb)
        np.testing.assert_allclose(wa@x,u,atol=1e-12)
        np.testing.assert_allclose(cosine_from_grams(wa,wb,x@x.T),cosine_matrix(u,v),atol=1e-12)
    def test_cross_gram(self):
        rng=np.random.default_rng(9);x=rng.normal(size=(12,6));y=rng.normal(size=(14,6))
        a={'a':{'g':[0,1]},'b':{'g':[2,3]}};c={'ntc':[8,9]}
        wa=response_weights(12,a,c);wb=response_weights(14,a,c)
        np.testing.assert_allclose(cosine_from_grams(wa,wb,x@x.T,y@y.T,x@y.T),cosine_matrix(wa@x,wb@y),atol=1e-12)

if __name__=='__main__':unittest.main()
