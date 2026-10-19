"""Tests for leakage, replicate weighting and denominator integrity in ML."""
import unittest

import numpy as np
import pandas as pd

from bio_visuals import morphology_ml as ml


def fixture():
    rng = np.random.RandomState(7)
    rows, meta = [], []
    for index in range(6):
        job, animal = str(100 + index), 'animal' + str(index)
        meta.append(dict(job=job, animal_id=animal,
                         running_group=ml.HOME if index < 3 else ml.RUNNING,
                         pla_arm='test arm', source_file=job + '.csv'))
        oid = 0
        for log_area, shape, count in [(-2, .15, 24 + index * 2), (-.7, .5, 30), (.5, .85, 24)]:
            for _ in range(count):
                oid += 1
                rows.append(dict(job=job, object_type='PLA', object_id=oid,
                                 area=10 ** rng.normal(log_area, .04),
                                 min_diameter=rng.normal(shape, .01), max_diameter=1.,
                                 chat_present=float(oid % (2 if index < 3 else 4) == 0),
                                 x=float(oid * 100), y=float(index * 500 + oid * 25)))
    return pd.DataFrame(rows), pd.DataFrame(meta)


class MorphologyMLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.objects, cls.metadata = fixture()
        cls.result = ml.run_analysis(cls.objects, cls.metadata, sample_per_animal=60, n_init=2)

    def test_outcome_condition_and_coordinates_do_not_enter_fit(self):
        changed = self.objects.copy()
        changed.chat_present = 1 - changed.chat_present
        changed.x = changed.x * 100 + 12345
        changed.y = -changed.y
        meta = self.metadata.copy()
        meta.running_group = meta.running_group.map({ml.HOME: ml.RUNNING, ml.RUNNING: ml.HOME})
        other = ml.run_analysis(changed, meta, sample_per_animal=60, n_init=2)
        self.assertEqual(self.result['definitions']['models'], other['definitions']['models'])
        pd.testing.assert_frame_equal(self.result['profiles'], other['profiles'])
        pd.testing.assert_series_equal(self.result['animals'].n, other['animals'].n)

    def test_equal_animal_fit_and_refit_scaler_excludes_heldout(self):
        for model in self.result['definitions']['models']:
            self.assertEqual(model['n_per_animal'], 60)
            self.assertEqual(len(model['training_animals']), 6)
        audit = self.result['fit_audit']
        for row in audit.itertuples():
            self.assertNotIn(row.omitted_animal, row.training_animals.split(','))
            self.assertEqual(row.n_training_animals, 5)
            self.assertEqual(row.training_n, 5 * row.n_per_animal)
        # Recompute one fold independently: scaler is fitted after omission.
        prepared, _ = ml._validate(self.objects, self.metadata)
        heldout = audit.iloc[0].omitted_animal
        train = prepared.loc[prepared.animal_id.ne(heldout)]
        sample, _ = ml._balanced_sample(train, 60, ml.SEED)
        self.assertAlmostEqual(audit.iloc[0].scaler_log_area_center, sample.log10_area.median())
        self.assertAlmostEqual(audit.iloc[0].scaler_shape_center, sample['shape'].median())

    def test_missing_shape_excluded_only_from_ml_and_area_sensitivity_retains_it(self):
        objects = self.objects.copy()
        objects.loc[0, ['min_diameter', 'max_diameter']] = 0.
        objects.loc[1, 'area'] = np.nan
        objects.loc[2, 'min_diameter'] = 0.  # Valid zero ratio, separately audited.
        original = objects.copy(deep=True)
        result = ml.run_analysis(objects, self.metadata, sample_per_animal=60, n_init=2)
        pd.testing.assert_frame_equal(objects, original)
        audit = result['audit'].set_index('job').loc['100']
        self.assertEqual(audit.n_excluded, 2)
        self.assertEqual(audit.both_zero_diameters_n, 1)
        self.assertEqual(audit.nonfinite_morphology_n, 1)
        self.assertEqual(audit.retained_zero_min_n, 1)
        self.assertEqual(sum(np.asarray(g['n_all']).sum() for g in result['spatial']), len(objects))
        self.assertEqual(sum(np.asarray(g['n_excluded']).sum() for g in result['spatial']), 2)
        area = result['area_only_standardized'].iloc[0]
        self.assertEqual(area.valid_pla_n, len(objects) - 1)
        self.assertEqual(area.morphology_pla_n, len(objects) - 2)
        self.assertLess(area.input_coverage_pct, 100)

    def test_common_mix_removes_pure_composition_difference(self):
        rows, raw = [], []
        for i in range(6):
            group = ml.HOME if i < 3 else ml.RUNNING
            counts = [10, 90] if i < 3 else [90, 10]
            for c, n, association in zip([1, 2], counts, [20, 80]):
                rows.append(dict(pla_arm='arm', animal_id=str(i), job=str(i),
                                 running_group=group, cluster=c, n=n,
                                 proportion_pct=n, association_pct=association))
            raw.append(dict(pla_arm='arm', animal_id=str(i), running_group=group,
                            raw_association_pct=74 if i < 3 else 26))
        standard, weights, individual = ml._standardize(pd.DataFrame(rows), pd.DataFrame(raw), 2, 10)
        self.assertAlmostEqual(standard[0]['raw_difference_pp'], 48)
        self.assertAlmostEqual(standard[0]['valid_difference_pp'], 48)
        self.assertAlmostEqual(standard[0]['standardized_difference_pp'], 0)
        self.assertEqual(standard[0]['retained_pct'], 100)
        self.assertTrue(all(abs(w['common_weight'] - .5) < 1e-12 for w in weights))
        # One unsupported cluster is omitted for ALL animals and reported.
        rows[0]['n'] = 9
        standard, weights, _ = ml._standardize(pd.DataFrame(rows), pd.DataFrame(raw), 2, 10)
        self.assertEqual(standard[0]['eligible_clusters'], 1)
        self.assertLess(standard[0]['retained_pct'], 100)
        self.assertFalse(weights[0]['eligible'])

    def test_missing_outcome_never_becomes_negative_or_vanishes(self):
        objects = self.objects.copy()
        objects.loc[0, 'chat_present'] = np.nan
        result = ml.run_analysis(objects, self.metadata, sample_per_animal=60, n_init=2)
        self.assertEqual(result['audit'].missing_association_n.sum(), 1)
        self.assertTrue(result['standardized'].raw_difference_pp.isna().all())
        missing = result['animals'].query('job == "100"').association_pct.isna().sum()
        self.assertEqual(missing, 1)
        self.assertLess(result['standardized'].retained_pct.max(), 100)

    def test_partitions_are_repeatable_and_ari_is_label_invariant(self):
        shuffled = self.objects.sample(frac=1, random_state=51).reset_index(drop=True)
        other = ml.run_analysis(shuffled, self.metadata, sample_per_animal=60, n_init=2)
        self.assertEqual(self.result['definitions']['models'], other['definitions']['models'])
        pd.testing.assert_frame_equal(self.result['animals'], other['animals'])
        self.assertEqual(ml._adjusted_rand([1, 1, 2, 2], [5, 5, 3, 3]), 1)
        self.assertAlmostEqual(ml._adjusted_rand([1, 1, 2, 2], [1, 2, 1, 2]), -.5)
        self.assertEqual(ml._adjusted_rand([1, 1, 1, 1], [1, 1, 2, 2]), 0)
        self.assertTrue(self.result['stability_clusters'].jaccard.between(0, 1).all())

    def test_native_spatial_percentages_conserve_valid_objects(self):
        valid_count = self.result['audit'].n_valid.sum()
        self.assertEqual(sum(np.asarray(grid['n']).sum() for grid in self.result['spatial']), valid_count)
        for grid in self.result['spatial']:
            counts = np.asarray(grid['n'])
            fractions = np.array([grid['cluster_fractions'][str(c)] for c in [1, 2, 3]], float)
            self.assertTrue(np.isnan(fractions[:, counts < 10]).all())
            np.testing.assert_allclose(fractions[:, counts >= 10].sum(axis=0), 100)
        self.assertEqual(len(self.result['review_targets']), 36)
        self.assertTrue(self.result['review_targets'].source_file.str.endswith('.csv').all())

    def test_excluded_qc_targets_preserve_provenance_and_do_not_change_fit(self):
        additions = []
        for job, count in [('100', 12), ('101', 8)]:
            source = self.objects.loc[self.objects.job.eq(job)].iloc[0].to_dict()
            for i in range(count):
                row = dict(source, object_id=1000 + i, min_diameter=0., max_diameter=0.,
                           x=float(i * 1100), y=float((i % 3) * 900))
                if i == count - 2:
                    row.update(area=0., min_diameter=.5, max_diameter=1.)
                elif i == count - 1:
                    row.update(min_diameter=-.2, max_diameter=1.)
                additions.append(row)
            additions.append(dict(source, object_id=2000, min_diameter=0., max_diameter=0., x=np.nan))
        objects = pd.concat([self.objects, pd.DataFrame(additions)], ignore_index=True)
        original = objects.copy(deep=True)
        result = ml.run_analysis(objects, self.metadata, sample_per_animal=60, n_init=2)
        pd.testing.assert_frame_equal(objects, original)
        targets = result['excluded_review_targets']
        self.assertEqual(targets.groupby('job').size().to_dict(), {'100': 6, '101': 6})
        self.assertFalse(targets[['job', 'object_id']].duplicated().any())
        self.assertTrue(np.isfinite(targets[['x', 'y', 'window_xmin', 'window_ymax']]).all().all())
        self.assertTrue(targets.section_excluded_missing_coordinate_n.eq(1).all())
        self.assertEqual(targets.groupby('job').section_excluded_n.first().to_dict(), {'100': 13, '101': 9})
        self.assertTrue(targets.reason.str.contains('not a prevalence estimate').all())
        for job, rows in targets.groupby('job'):
            self.assertEqual(set(rows.exclusion_reason), {'Both reported diameters are zero',
                                                         'Nonpositive PLA area', 'Negative minimum diameter'})
            self.assertEqual(rows[['native_bin_x', 'native_bin_y']].drop_duplicates().shape[0], 6)
            self.assertTrue(rows.source_file.eq(job + '.csv').all())
        source = objects.set_index(['job', 'object_id'])
        for row in targets.itertuples():
            original_row = source.loc[(row.job, row.object_id)]
            self.assertEqual((row.x, row.y), (original_row.x, original_row.y))
            self.assertEqual((row.min_diameter, row.max_diameter),
                             (original_row.min_diameter, original_row.max_diameter))
            self.assertEqual(row.window_xmin, row.x - 500)
            self.assertEqual(row.window_ymax, row.y + 500)
        # Adding objects ineligible for morphology must not alter its fit or the
        # original eligible-object image-review queue.
        only_morphology = lambda r: [m for m in r['definitions']['models'] if m['feature_set'] == 'area_and_shape']
        self.assertEqual(only_morphology(self.result), only_morphology(result))
        pd.testing.assert_frame_equal(self.result['review_targets'], result['review_targets'])
        prepared, _ = ml._validate(objects.sample(frac=1, random_state=31), self.metadata)
        pd.testing.assert_frame_equal(targets, ml._excluded_review_targets(prepared))
        prepared.chat_present = 1 - prepared.chat_present
        unchanged = ml._excluded_review_targets(prepared)
        pd.testing.assert_frame_equal(targets.drop(columns='chat_present'), unchanged.drop(columns='chat_present'))

    def test_excluded_qc_empty_and_unlocatable_records_are_not_invented(self):
        empty = self.result['excluded_review_targets']
        self.assertTrue(empty.empty)
        self.assertEqual(list(empty.columns), ml.EXCLUDED_REVIEW_COLUMNS)
        changed = self.objects.copy()
        changed.loc[0, ['min_diameter', 'max_diameter']] = 0.
        changed.loc[0, 'x'] = np.nan
        prepared, _ = ml._validate(changed, self.metadata)
        targets = ml._excluded_review_targets(prepared)
        self.assertTrue(targets.empty)
        self.assertEqual(list(targets.columns), ml.EXCLUDED_REVIEW_COLUMNS)


if __name__ == '__main__':
    unittest.main()
