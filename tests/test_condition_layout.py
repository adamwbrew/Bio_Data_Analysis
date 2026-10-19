"""Scientific safeguards for the condition-versus-image-layout audit."""
import unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bio_visuals.condition_layout import layout_audit, layout_figure
from bio_visuals.handoff import ARMS


class ConditionLayoutTests(unittest.TestCase):
    def fixture(self):
        metadata, objects = [], []
        for arm in ARMS.values():
            for animal, group, x, flags in [('21F', 'Home cage', 100., [True, True]),
                                             ('22F', 'Running', 0., [True, False])]:
                job = arm + animal
                metadata.append(dict(job=job, animal_id=animal, running_group=group, pla_arm=arm,
                                     image_id='source-one', image_label='Image 1'))
                for i, flag in enumerate(flags):
                    objects.append(dict(job=job, xmin=x+i, xmax=x+i+4, ymin=10., ymax=14.,
                                        x=x+i+2, y=12., object_type='PLA', chat_present=flag))
                # Broad ChAT object tests that bounds include ALL classes, not PLA only.
                objects.append(dict(job=job, xmin=x-2, xmax=x+8, ymin=8., ymax=18.,
                                    x=x+3, y=13., object_type='ChAT', chat_present=True))
        return pd.DataFrame(objects), pd.DataFrame(metadata)

    def test_extents_use_original_boxes_and_contrasts_have_requested_sign(self):
        objects, metadata = self.fixture()
        extent, contrasts, summary = layout_audit(objects, metadata)
        home = extent[extent.running_group.eq('Home cage')]
        np.testing.assert_array_equal(home[['x_min','x_max','y_min','y_max']], [[98.,108.,8.,18.]]*2)
        np.testing.assert_array_equal(contrasts.home_minus_running, [50., 50.])
        self.assertEqual(summary['condition_x_order'], 'Running lower-X; Home cage higher-X')
        self.assertTrue(summary['complete_condition_x_separation'])
        self.assertEqual(summary['positive_within_image_contrasts'], 2)
        self.assertEqual(len(extent), 4)
        self.assertTrue(contrasts.status.str.contains('not a matched experimental pair').all())

    def test_overlap_is_not_called_separation_and_opposite_effect_is_preserved(self):
        objects, metadata = self.fixture()
        for job in metadata.loc[metadata.running_group.eq('Home cage'), 'job']:
            use = objects.job.eq(job)
            objects.loc[use, ['xmin','xmax']] -= 98
            objects.loc[use & objects.object_type.eq('PLA'), 'chat_present'] = False
        _, contrasts, summary = layout_audit(objects, metadata)
        self.assertFalse(summary['complete_condition_x_separation'])
        np.testing.assert_array_equal(contrasts.home_minus_running, [-50., -50.])
        self.assertEqual(summary['negative_within_image_contrasts'], 2)

    def test_conflicting_assignment_duplicates_and_unknown_jobs_rejected(self):
        objects, metadata = self.fixture()
        bad = metadata.copy(); bad.loc[0, 'running_group'] = 'Running'
        with self.assertRaisesRegex(ValueError, 'Condition assignments'):
            layout_audit(objects, bad)
        with self.assertRaisesRegex(ValueError, 'unique section'):
            layout_audit(objects, pd.concat([metadata, metadata.iloc[[0]]]))
        bad = objects.copy(); bad.loc[0, 'job'] = 'unreviewed-duplicate'
        with self.assertRaisesRegex(ValueError, 'primary section metadata'):
            layout_audit(bad, metadata)

    def test_missing_association_stays_unknown_and_invalid_coordinates_rejected(self):
        objects, metadata = self.fixture()
        objects['chat_present'] = objects.chat_present.astype('object')
        objects.loc[0, 'chat_present'] = np.nan
        _, contrasts, summary = layout_audit(objects, metadata)
        self.assertTrue(np.isnan(contrasts.iloc[0].home_minus_running))
        self.assertEqual(summary['finite_within_image_contrasts'], 1)
        bad = objects.copy(); bad.loc[0, 'xmin'] = np.inf
        with self.assertRaisesRegex(ValueError, 'finite numeric'):
            layout_audit(bad, metadata)
        bad = objects.copy(); bad.loc[0, 'xmax'] = -1
        with self.assertRaisesRegex(ValueError, 'lower bounds'):
            layout_audit(bad, metadata)

    def test_missing_condition_within_image_is_not_called_a_pair(self):
        objects, metadata = self.fixture()
        metadata.loc[0, 'image_id'] = 'other-source'
        metadata.loc[0, 'image_label'] = 'Image 2'
        _, contrasts, summary = layout_audit(objects, metadata)
        self.assertFalse(summary['complete_condition_x_separation'])
        self.assertEqual(contrasts.home_minus_running.isna().sum(), 2)
        self.assertEqual(summary['finite_within_image_contrasts'], 1)

    def test_figure_keeps_arms_separate_and_draws_both_condition_extents(self):
        objects, metadata = self.fixture()
        extent, _, _ = layout_audit(objects, metadata)
        fig = layout_figure(extent)
        self.assertEqual(len(fig.axes), 2)
        for axis in fig.axes:
            self.assertEqual(len(axis.lines), 2)
            self.assertIn('TDP43', axis.get_title(loc='left'))
        plt.close(fig)


if __name__ == '__main__':
    unittest.main()
