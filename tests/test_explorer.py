"""Measurement and alignment boundaries for the notebook explorer."""
import copy
import json
import unittest
import numpy as np
import pandas as pd
from bio_visuals import handoff as h
from bio_visuals.explorer import explorer_payload, _grid_payload


class ExplorerTests(unittest.TestCase):
    def fixture(self):
        metadata=[];reviews=[];grids={}
        for arm in h.ARMS.values():
            for i,animal in enumerate(h.ANIMAL_ORDER):
                job=arm+animal
                metadata.append(dict(job=job,animal_id=animal,running_group=h.GROUPS[animal],
                                     pla_arm=arm,source_file=job+'.csv'))
                reviews.append(dict(job=job,animal_id=animal,running_group=h.GROUPS[animal],pla_arm=arm,
                                    landmark_reviewed=False,orientation_reviewed=False,
                                    scale_common_confirmed=False,support_reviewed=False,orientation='horizontal',
                                    candidate=dict(centroid_native=[7,12],polygon_native=[[0,0],[14,0],[7,36]]),
                                    reviewed_polygon_native=None,reviewer='',second_reviewer='',support_qc='uncertain'))
                f=pd.DataFrame(dict(x=[1,2,12],y=[1,2,1],xmin=[0,1,11],xmax=[2,3,13],
                                    ymin=[0,1,0],ymax=[2,3,2],object_type=['PLA','PLA','ChAT'],
                                    area=[2.,5.,3.],overlap_area=[.5,0,0],chat_present=[True,False,True]))
                grids[job]=h.build_grid(f,10,10)
        return grids,pd.DataFrame(metadata),reviews

    def test_payload_preserves_native_coordinates_counts_and_area(self):
        grids,metadata,reviews=self.fixture()
        result=explorer_payload(grids,grids,grids,grids,metadata,reviews)
        json.dumps(result,allow_nan=False)
        for s in result['sections']:
            g=s['grids']['10']
            self.assertEqual(g['x'],[5.,15.])
            self.assertEqual(np.sum(g['n_pla']),2)
            self.assertEqual(np.sum(g['n_chat']),1)
            self.assertEqual(np.sum(g['n_positive']),1)
            self.assertEqual(np.sum(g['pla_area']),7)
            self.assertEqual(np.sum(g['overlap_area']),.5)

    def test_unapproved_reviews_cannot_supply_reviewed_transforms(self):
        grids,metadata,reviews=self.fixture()
        result=explorer_payload(grids,grids,grids,grids,metadata,reviews)
        for s in result['sections']:
            self.assertFalse(s['ready'])
            self.assertIsNone(s['reviewed_center'])
            self.assertIsNone(s['reviewed_flip'])
            self.assertEqual(s['candidate_center'],[7,12])

    def test_reviewed_transform_requires_all_flags_and_second_flip_review(self):
        grids,metadata,reviews=self.fixture()
        r=reviews[0]
        r.update(landmark_reviewed=True,orientation_reviewed=True,scale_common_confirmed=True,
                 support_reviewed=True,support_qc='complete',reviewer='A',second_reviewer='B',
                 reviewed_polygon_native=[[0,0],[4,0],[4,6],[0,6]])
        result=explorer_payload(grids,grids,grids,grids,metadata,reviews)
        self.assertTrue(result['sections'][0]['ready'])
        self.assertEqual(result['sections'][0]['reviewed_center'],[2,3])
        self.assertEqual(result['sections'][0]['reviewed_flip'],'horizontal')
        r['second_reviewer']=''
        self.assertFalse(explorer_payload(grids,grids,grids,grids,metadata,reviews)['sections'][0]['ready'])

    def test_display_order_restarts_per_arm_and_missing_candidates_stay_missing(self):
        grids,metadata,reviews=self.fixture();reviews[0]['candidate']=None
        sections=explorer_payload(grids,grids,grids,grids,metadata,reviews)['sections']
        self.assertIsNone(sections[0]['candidate_center'])
        for arm in h.ARMS.values():
            a=[s for s in sections if s['arm']==arm]
            self.assertEqual([s['animal'] for s in a],h.ANIMAL_ORDER)
            self.assertEqual([s['order'] for s in a],list(range(6)))
            self.assertTrue(all('depth' not in s for s in a))

    def test_payload_rejects_missing_or_negative_measurements(self):
        grids,_,_=self.fixture();g=copy.deepcopy(next(iter(grids.values())))
        g['n_pla'][0,0]=np.nan
        with self.assertRaises(ValueError):_grid_payload(g)


if __name__=='__main__':unittest.main()
