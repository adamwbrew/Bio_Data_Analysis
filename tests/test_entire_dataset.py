"""Scientific invariants for expanded data, duplicate handling and visual support."""
import csv
from pathlib import Path
import tempfile
import unittest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bio_visuals.data import SOURCE_COLUMNS
from bio_visuals.entire_dataset import (discover_exports,load_complete_dataset,load_export,
    signature_overlap,run_summary,type_summary,metadata_template,biological_summaries)
from bio_visuals.complete_plots import (native_grids,grid_table,section_figure,
    shared_count_scales,aligned_rerun_grids,overview_figure,layout_figure)

IMAGE='20260515_173745_1_ZAbMVm_Ellie_PLA_05152026_1'


def write_export(root,job,region='22F_TDP43-Pfkp_tight',rows=4):
    path=Path(root)/f'{IMAGE}_backSUB.ome.ome.tiff_job{job}PLA and ChaT Object Colocalization FL v3.0.0_object_results.csv'
    with path.open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.writer(handle);writer.writerow(SOURCE_COLUMNS)
        for i in range(rows):
            typ='PLA - TRITC_E1' if i%2==0 else 'ChAt-50 - FITC_E1'
            writer.writerow([f'H:\\{IMAGE}\\{IMAGE}_backSUB.ome.ome.tiff',region,'algorithm',i,typ,
                             1,100,1,0,1,.2,0,1,.5,1,1,10+i*3,11+i*3,10,11])
    return path


class CompleteDatasetTests(unittest.TestCase):
    def tearDown(self):plt.close('all')

    def test_only_explicit_replacement_is_selected_and_excluded_export_is_audited(self):
        with tempfile.TemporaryDirectory() as root:
            write_export(root,'4978',rows=2);write_export(root,'5403',rows=4)
            objects,manifest=load_complete_dataset(root)
            self.assertEqual(set(objects.job),{'5403'})
            self.assertEqual(manifest.set_index('job').loc['4978','duplicate_status'],'interrupted_superseded')
            self.assertIn('overlap_presence_mismatch',manifest.columns)
            self.assertTrue(manifest.sha256.str.len().eq(64).all())
            short=load_export(next(Path(root).glob('*job4978*')))
            audit=signature_overlap(short,objects)
            self.assertEqual(audit['short_unmatched_rows'],0)
            self.assertEqual(audit['full_additional_rows'],2)
            # IDs do not establish object matches; changed measurements do.
            changed=short.copy();changed['object_id']+=50
            self.assertEqual(signature_overlap(changed,objects)['identical_measurement_rows'],2)
            changed.loc[0,'area']+=.1
            self.assertEqual(signature_overlap(changed,objects)['identical_measurement_rows'],1)

    def test_unreviewed_duplicate_and_missing_replacement_fail(self):
        with tempfile.TemporaryDirectory() as root:
            write_export(root,'4978')
            with self.assertRaisesRegex(ValueError,'requires replacement'):discover_exports(root)
        with tempfile.TemporaryDirectory() as root:
            write_export(root,'100');write_export(root,'101',rows=6)
            with self.assertRaisesRegex(ValueError,'Unreviewed repeated'):discover_exports(root)

    def test_bad_flags_ids_and_later_metadata_are_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            path=write_export(root,'5403')
            base=pd.read_csv(path,encoding='utf-8-sig')
            for column,value,pattern in [('PLA - TRITC_E1 present',2,'Nonbinary'),
                    ('Object Id',.5,'Object Id'),('Analysis Region','21F_TDP43-Pfkp_tight','Inconsistent')]:
                edited=base.copy()
                if column == 'Object Id':
                    edited[column]=edited[column].astype(float)
                edited.loc[1,column]=value
                edited.to_csv(path,index=False,encoding='utf-8-sig')
                with self.assertRaisesRegex(ValueError,pattern):load_export(path)

    def test_native_grids_conserve_counts_mask_low_support_and_preserve_zero_fraction(self):
        with tempfile.TemporaryDirectory() as root:
            path=write_export(root,'5403',rows=6);o=load_export(path)
            o.loc[o.object_type=='PLA','chat_present']=False
            grids=native_grids(o,bin_width=3,min_support=1);g=grids['5403']
            self.assertEqual(g['pla_count'].sum(),3);self.assertEqual(g['chat_count'].sum(),3)
            self.assertTrue(np.isnan(g['fraction'][g['pla_count']==0]).all())
            np.testing.assert_array_equal(g['fraction'][g['pla_count']>0],0)
            masked=native_grids(o,bin_width=3,min_support=2)['5403']
            self.assertTrue(np.isnan(masked['fraction']).all())
            table=grid_table(grids);self.assertEqual(table.pla_count.sum(),3)
            fig=section_figure(g,g,shared_count_scales(grids));self.assertEqual(len(fig.axes),8)
            for width in [0,-1,np.nan]:
                with self.assertRaises(ValueError):native_grids(o,bin_width=width)

    def test_duplicate_multiset_matching_respects_multiplicity(self):
        with tempfile.TemporaryDirectory() as root:
            o=load_export(write_export(root,'5403'))
            short=pd.concat([o.iloc[[0]],o.iloc[[0]]],ignore_index=True)
            self.assertEqual(signature_overlap(short,o)['identical_measurement_rows'],1)

    def test_equal_export_mean_is_distinct_from_detection_weighting(self):
        s=pd.DataFrame({'pla_type':['Hk1','Hk1'],'pla_objects':[10,90],'pla_with_chat':[10,0],
                        'pla_association_pct':[100,0],'area_weighted_overlap_pct':[50,10]})
        result=type_summary(s).iloc[0]
        self.assertEqual(result.mean_export_association_pct,50)
        self.assertEqual(result.pooled_object_association_pct,10)

    def test_group_means_require_metadata_and_weight_animals_equally(self):
        s=pd.DataFrame(dict(job=['1','2','3'],condition=['A','A','B'],pla_type=['Hk1']*3,
            image_id=['image']*3,pla_association_pct=[0,100,100],area_weighted_overlap_pct=[0,100,100]))
        metadata=metadata_template(s)
        animals,groups,missing=biological_summaries(s,metadata)
        self.assertTrue(animals.empty);self.assertEqual(len(missing),4)
        metadata['animal_id']=['a','a','b'];metadata['running_group']='running'
        metadata['section_id']=['s1','s2','s3'];metadata['anatomical_region']='region'
        animals,groups,missing=biological_summaries(s,metadata)
        self.assertEqual(len(animals),2);self.assertEqual(groups.mean_association_pct.iloc[0],75)
        metadata.loc[0,'running_group']='control'
        with self.assertRaisesRegex(ValueError,'longitudinal'):biological_summaries(s,metadata)

    def test_summary_and_native_layout_show_exports_without_invented_registration(self):
        with tempfile.TemporaryDirectory() as root:
            write_export(root,'5403');write_export(root,'4997',region='22F_Hk1-TDP43_tight')
            objects,manifest=load_complete_dataset(root);s=run_summary(objects)
            self.assertEqual(len(s),2);self.assertEqual(s.all_objects.sum(),8)
            self.assertEqual(len(overview_figure(s).axes),3)
            self.assertEqual(len(layout_figure(s).axes),1)


if __name__=='__main__':unittest.main()
