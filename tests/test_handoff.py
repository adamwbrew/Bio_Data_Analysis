"""Scientific invariants for the raw-data spatial handoff (small synthetic data)."""
import copy
import unittest
import numpy as np
import pandas as pd
from bio_visuals import handoff as h


class SpatialHandoffTests(unittest.TestCase):
    def test_translation_and_flips_preserve_distances_and_center(self):
        xy=np.array([[2.,3.],[5.,7.],[-1.,8.]])
        center=np.array([4.,6.])
        for flip in ('none','horizontal','vertical','both'):
            oriented,aligned=h.transform_xy(xy,center,flip)
            np.testing.assert_allclose(np.linalg.norm(np.diff(aligned,axis=0),axis=1),np.linalg.norm(np.diff(xy,axis=0),axis=1))
            np.testing.assert_allclose(h.transform_xy(center,center,flip)[1],[0,0])
        with self.assertRaises(ValueError):h.transform_xy(xy,center,'rotate 90')

    def test_polygon_centroid_is_area_weighted_not_vertex_mean(self):
        # Same rectangle with extra vertices on one side: centroid stays (1,1).
        poly=[[0,0],[2,0],[2,.1],[2,.2],[2,2],[0,2]]
        np.testing.assert_allclose(h.polygon_centroid(poly),[1,1])
        with self.assertRaises(ValueError):h.polygon_centroid([[0,0],[1,0],[2,0]])

    def test_raw_grid_keeps_zeros_distinct_from_unmeasured(self):
        f=pd.DataFrame(dict(x=[1.,2.,12.],y=[1.,2.,1.],xmin=[0,1,11],xmax=[2,3,13],
                            ymin=[0,1,0],ymax=[2,3,2],object_type=['PLA','PLA','ChAT'],
                            chat_present=[False,False,True],area=[1.,2.,3.],overlap_area=[0.,0.,0.]))
        g=h.build_grid(f,width=10,min_pla=2)
        self.assertEqual(g['n_all'].sum(),3)
        self.assertEqual(g['n_pla'].sum(),2)
        self.assertEqual(g['fraction'][0,0],0)  # two measured negative PLA
        self.assertTrue(np.isnan(g['fraction'][0,1]))  # ChAT only; no PLA denominator

    def test_positive_object_area_is_not_actual_overlap_area(self):
        f=pd.DataFrame(dict(x=[1.,2.],y=[1.,2.],xmin=[0,1],xmax=[2,3],ymin=[0,1],ymax=[2,3],
                            object_type=['PLA','PLA'],chat_present=[True,False],area=[10.,20.],overlap_area=[2.,0.]))
        g=h.build_grid(f,width=10,min_pla=1)
        self.assertEqual(g['positive_pla_area'].sum(),10)
        self.assertEqual(g['overlap_area'].sum(),2)
        self.assertEqual(g['pla_area'].sum(),30)

    def test_equal_animal_mean_and_missing_coverage(self):
        values=np.array([[[10.,np.nan]],[[90.,np.nan]],[[50.,20.]]])
        mean,n=h.equal_animal_maps(values,minimum=3)
        self.assertEqual(mean[0,0],50.)
        self.assertTrue(np.isnan(mean[0,1]))
        np.testing.assert_array_equal(n,[[3,1]])

    def fixture(self):
        records=[];rows=[]
        for arm in h.ARMS.values():
            for i,(animal,group) in enumerate(h.GROUPS.items()):
                job=arm+animal
                records.append(dict(job=job,animal_id=animal,pla_arm=arm,running_group=group,
                                    reviewed_polygon_native=[[-1,-1],[1,-1],[1,1],[-1,1]],
                                    landmark_reviewed=True,orientation='none',orientation_reviewed=True,
                                    scale_common_confirmed=True,support_reviewed=True,support_qc='complete',
                                    reviewer='Synthetic test',second_reviewer=''))
                # Both arms deliberately have opposite effects, with unequal raw counts.
                positive=(group=='Running') == (arm=='TDP43–HK1')
                for k in range(i+1):rows.append(dict(job=job,x=10.,y=10.,object_type='PLA',chat_present=positive))
        return records,pd.DataFrame(rows)

    def test_review_gate_prevents_accidental_group_maps(self):
        records,objects=self.fixture()
        for field in ('landmark_reviewed','orientation_reviewed','scale_common_confirmed','support_reviewed'):
            bad=copy.deepcopy(records);bad[0][field]=False
            with self.assertRaisesRegex(ValueError,'withheld'):h.registered_fraction_maps(objects,bad,min_pla=1)
        bad=copy.deepcopy(records);bad[0]['orientation']='horizontal'
        with self.assertRaisesRegex(ValueError,'withheld'):h.registered_fraction_maps(objects,bad,min_pla=1)
        bad=copy.deepcopy(records);bad[0]['landmark_reviewed']='false'
        with self.assertRaisesRegex(ValueError,'JSON boolean'):h.review_status(bad)

    def test_surface_selection_preserves_each_sections_xy_ratio(self):
        from bio_visuals.handoff_plots import interactive_surface
        records,_=self.fixture();metadata=pd.DataFrame(records)
        grids={}
        for i,r in enumerate(records):
            grids[r['job']]=dict(x_edges=np.arange(4)*10.,y_edges=np.arange(3)*(i+1)*10.,
                                 n_pla=np.ones((2,3)),n_chat=np.ones((2,3)),
                                 fraction=np.full((2,3),np.nan),support=np.ones((2,3),bool))
        fig=interactive_surface(grids,metadata,'TDP43–HK1')
        jobs=metadata[metadata.pla_arm=='TDP43–HK1'].set_index('animal_id')
        for step in fig.layout.sliders[0].steps:
            g=grids[jobs.loc[step.label,'job']]
            expected=np.ptp(g['y_edges'])/np.ptp(g['x_edges'])
            self.assertAlmostEqual(step.args[1]['scene.aspectratio']['y'],expected)

    def test_reviewed_maps_never_pool_the_two_arms(self):
        records,objects=self.fixture()
        result=h.registered_fraction_maps(objects,records,min_pla=1)
        hk=result['arms']['TDP43–HK1']['difference'];pf=result['arms']['TDP43–PFKP']['difference']
        np.testing.assert_allclose(hk[np.isfinite(hk)],[100.])
        np.testing.assert_allclose(pf[np.isfinite(pf)],[-100.])
        self.assertEqual(np.isfinite(hk).sum(),1)

    def test_condition_summary_rejects_repeated_animal(self):
        records,_=self.fixture()
        s=pd.DataFrame(records);s['pla_association_pct']=np.arange(len(s))
        with self.assertRaises(ValueError):h.group_contrasts(pd.concat([s,s.iloc[[0]]]))

    def test_notebook_python_compiles(self):
        from pathlib import Path
        import nbformat
        nb=nbformat.read(Path(__file__).resolve().parents[1]/'PLA_ChAT_spatial_handoff_analysis.ipynb',as_version=4)
        for i,c in enumerate(nb.cells):
            if c.cell_type=='code':compile(c.source,f'notebook cell {i}','exec')

if __name__=='__main__':unittest.main()
