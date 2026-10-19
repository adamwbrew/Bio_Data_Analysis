import copy
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from bio_visuals import condition_comparison as cc, feature_alignment as fa, handoff as h
from bio_visuals.explorer import explorer_payload
from tests.test_condition_comparison import fixture
from tests import test_explorer

class NormalizationTests(unittest.TestCase):
    def test_two_directions_reverse_estimate_and_interval_not_relative_baseline(self):
        o,m=fixture([20,40,60,0,10,20]);r=cc.comparison_results(o,m)
        for _,d in r['directions'].groupby('pla_arm'):
            a,b=d.iloc[0],d.iloc[1]
            self.assertAlmostEqual(a.difference_pp,-b.difference_pp)
            self.assertAlmostEqual(a.leave_pair_out_min,-b.leave_pair_out_max)
            self.assertAlmostEqual(a.leave_pair_out_max,-b.leave_pair_out_min)
            self.assertAlmostEqual(a.ratio*b.ratio,1)
            self.assertNotAlmostEqual(a.relative_change_pct,-b.relative_change_pct)
            self.assertAlmostEqual(a.between_animal_se_pp,np.sqrt(400/3+100/3))
    def test_denominators_expose_unequal_pooling_and_keep_equal_animals(self):
        o,m=fixture([90,20,20,10,10,10],[10,100,100,100,100,100]);r=cc.comparison_results(o,m)
        d=r['denominators'];self.assertTrue(np.allclose(d.equal_animal_weight_pct,100/3))
        self.assertTrue(np.allclose(d.groupby(['pla_arm','running_group']).pooled_object_weight_pct.sum(),100))
        p=r['normalization'].iloc[:4]
        self.assertAlmostEqual(p.iloc[0].home_mean,130/3)
        self.assertAlmostEqual(p.iloc[3].home_mean,100*49/210)
    def test_missing_association_and_overlap_do_not_disappear_in_audit(self):
        o,m=fixture();o['chat_present']=o.chat_present.astype(float);o.loc[0,'chat_present']=np.nan;o.loc[0,'overlap_area']=np.nan
        r=cc.comparison_results(o,m);n=r['normalization'].iloc[:4]
        self.assertTrue(np.isnan(n.iloc[0].home_mean));self.assertTrue(np.isnan(n.iloc[2].home_mean));self.assertTrue(np.isnan(n.iloc[3].home_mean))
        json.dumps(cc.comparison_payload(r),allow_nan=False)
    def test_within_animal_arm_difference_uses_same_identity(self):
        o,m=fixture(inverse_second_arm=True);r=cc.comparison_results(o,m)
        p=r['paired_arms'].set_index('animal_id')
        self.assertEqual(p.loc['21F','pfkp_minus_hk1_pp'],-60)
        self.assertEqual(p.loc['22F','pfkp_minus_hk1_pp'],100)

class AlignmentTests(unittest.TestCase):
    def test_translation_and_flips_preserve_all_pairwise_distances(self):
        x=np.array([[1.,2.],[4.,7.],[3.,-2.]])
        for flip in fa.FLIPS:
            a=h.transform_xy(x,[2,4],flip)[1]
            np.testing.assert_allclose(np.linalg.norm(a[:,None]-a,axis=2),np.linalg.norm(x[:,None]-x,axis=2))
    def test_boundary_score_is_symmetric_and_exact_matches_are_zero(self):
        a=np.array([[0.,0.],[1,0],[0,2],[3,4]])
        self.assertEqual(fa.boundary_score(a,a),(0.,0.))
        self.assertEqual(fa.boundary_score(a,a+[1,2]),fa.boundary_score(a+[1,2],a))
    def test_partial_boundary_trim_tolerates_small_missing_piece_but_reports_tail(self):
        a=np.column_stack([np.arange(100.),np.zeros(100)])
        b=np.vstack([a,np.array([[1000.,0.]]*15)])
        score,tail=fa.boundary_score(a,b)
        self.assertEqual(score,0);self.assertGreater(tail,800)
    def test_proposals_never_unlock_reviewed_coordinates(self):
        grids,m,reviews=test_explorer.ExplorerTests().fixture();original=copy.deepcopy(reviews)
        proposals=[dict(job=r['job'],center=[0.,0.],flip='vertical',approved=False) for r in reviews]
        p=explorer_payload(grids,grids,grids,grids,m,reviews,proposals=proposals)
        self.assertTrue(all(s['proposed'] for s in p['sections']))
        self.assertTrue(all(not s['ready'] and s['reviewed_center'] is None for s in p['sections']))
        self.assertEqual(reviews,original)
    def test_annotation_hash_mismatch_fails_closed(self):
        grids,m,reviews=test_explorer.ExplorerTests().fixture();m['sha256']='expected'
        annotations=[dict(job=r['job'],animal_id=r['animal_id'],pla_arm=r['pla_arm'],source_sha256='wrong',polygon_native=[[0,0],[4,0],[4,4],[0,4]]) for r in reviews]
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'annotations.json';path.write_text(json.dumps(annotations))
            with self.assertRaisesRegex(ValueError,'hash mismatch'):fa.alignment_proposals(grids,m,path)


class OrientationCueTests(unittest.TestCase):
    def test_internal_cue_selects_vertical_family_without_moving_anchor(self):
        self.assertEqual(fa.allowed_flips([2,3],[10,20],[0,8],[9,15]),('vertical','both'))
        self.assertEqual(fa.allowed_flips([2,3],[10,20],[0,1],[9,15]),('none','horizontal'))
    def test_missing_or_horizontal_cue_preserves_all_candidates(self):
        self.assertEqual(fa.allowed_flips([0,0],[0,0],None,[1,2]),fa.FLIPS)
        self.assertEqual(fa.allowed_flips([0,0],[0,0],[1,0],[1,2]),fa.FLIPS)

if __name__=='__main__':unittest.main()
