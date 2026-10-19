"""Scientific tests for condition subtraction and size-composition sensitivity."""
import json
import unittest
from unittest import mock

import numpy as np
import pandas as pd

from bio_visuals import condition_comparison as cc
from bio_visuals.handoff import ANIMAL_ORDER, ARMS, GROUPS


def fixture(rates=None, counts=None, inverse_second_arm=False):
    rates = rates or [80, 90, 100, 0, 10, 20]
    counts = counts or [100] * 6
    rows, metadata = [], []
    for arm_index, arm in enumerate(ARMS.values()):
        for index, animal in enumerate(ANIMAL_ORDER):
            job = f'{arm_index}_{animal}'
            metadata.append(dict(job=job, animal_id=animal, pla_arm=arm, running_group=GROUPS[animal]))
            count = counts[index]
            rate = 100 - rates[index] if inverse_second_arm and arm_index else rates[index]
            positive = round(rate * count / 100)
            for k in range(count):
                rows.append(dict(job=job, object_type='PLA', area=.2, chat_present=k < positive,
                                 overlap_area=.02 if k < positive else 0.))
            rows.append(dict(job=job, object_type='ChAT', area=1., chat_present=True, overlap_area=0.))
    return pd.DataFrame(rows), pd.DataFrame(metadata)


class ConditionComparisonTests(unittest.TestCase):
    def test_equal_animal_means_direction_and_arm_separation(self):
        objects, metadata = fixture([90, 20, 20, 10, 10, 10], [10, 100, 100, 100, 100, 100], True)
        result = cc.comparison_results(objects, metadata)
        effects = result['contrasts'].query('metric == "association_pct"').set_index('pla_arm')
        first, second = list(ARMS.values())
        self.assertAlmostEqual(effects.loc[first, 'home_mean'], 130 / 3)
        self.assertAlmostEqual(effects.loc[first, 'home_minus_running'], 100 / 3)
        self.assertAlmostEqual(effects.loc[second, 'home_minus_running'], -100 / 3)
        pooled = 100 * (9 + 20 + 20) / 210
        self.assertNotAlmostEqual(effects.loc[first, 'home_mean'], pooled)
        self.assertEqual(set(effects.n_home), {3})
        self.assertEqual(set(effects.n_running), {3})

    def test_nine_delete_pair_sensitivities_are_not_independent_pairs(self):
        objects, metadata = fixture([20, 40, 60, 0, 10, 20])
        result = cc.comparison_results(objects, metadata)
        arm = next(iter(ARMS.values()))
        effect = result['contrasts'].query('pla_arm == @arm and metric == "association_pct"').iloc[0]
        self.assertEqual(effect.home_minus_running, 30)
        self.assertEqual(effect.leave_pair_out_min, 15)
        self.assertEqual(effect.leave_pair_out_max, 45)
        cross = result['cross_animal'].query('pla_arm == @arm')
        self.assertEqual(len(cross), 9)
        self.assertEqual(cross.home_animal.nunique(), 3)
        self.assertEqual(cross.running_animal.nunique(), 3)
        self.assertEqual(cross.home_minus_running.min(), 0)
        self.assertEqual(cross.home_minus_running.max(), 60)
        self.assertIn('not nine independent', result['definitions']['cross_animal'])

    def test_exact_enumeration_only_primary_metric_and_minimum_two_sided_p(self):
        objects, metadata = fixture()
        result = cc.comparison_results(objects, metadata)
        for arm in ARMS.values():
            null = result['exact_allocations'].query('pla_arm == @arm')
            self.assertEqual(len(null), 20)
            self.assertEqual(null.is_observed.sum(), 1)
            self.assertEqual(len(set(null.home_animals)), 20)
            np.testing.assert_allclose(sorted(null.home_minus_running), sorted(-null.home_minus_running))
        np.testing.assert_allclose(result['exact_summary'].two_sided_p, [.1, .1])
        self.assertIn('exchangeable', result['definitions']['exact_test'])
        self.assertEqual(result['definitions']['primary_metric'], 'association_pct')

    def test_fixed_common_size_mix_distinguishes_composition_from_within_size_rates(self):
        rows, metadata = [], []
        for arm_index, arm in enumerate(ARMS.values()):
            for animal in ANIMAL_ORDER:
                job = f'{arm_index}_{animal}'
                group = GROUPS[animal]
                metadata.append(dict(job=job, animal_id=animal, running_group=group, pla_arm=arm))
                counts = [10, 90] if group == cc.HOME else [90, 10]
                for area, count, fraction in zip([.02, 1.5], counts, [.2, .8]):
                    for k in range(count):
                        rows.append(dict(job=job, object_type='PLA', area=area, chat_present=k < round(count * fraction), overlap_area=0.))
        result = cc.comparison_results(pd.DataFrame(rows), pd.DataFrame(metadata))
        summary = result['standardization']
        np.testing.assert_allclose(summary.query('estimator == "raw"').home_minus_running, [48, 48])
        np.testing.assert_allclose(summary.query('estimator == "restricted"').home_minus_running, [48, 48])
        np.testing.assert_allclose(summary.query('estimator == "standardized"').home_minus_running, [0, 0], atol=1e-12)
        selected = result['standard_weights'].query('included')
        np.testing.assert_allclose(selected.common_weight, [.5, .5, .5, .5])
        self.assertTrue(summary.min_retained_pct.eq(100).all())
        self.assertEqual(len(result['standard_weights'].query('not included')), 6)

    def test_support_rule_requires_every_animal_and_restricted_retention_is_reported(self):
        objects, metadata = fixture()
        first_job = objects.job.iloc[0]
        # Nine observations in the small stratum means this stratum is excluded
        # globally, although all its measurements are valid.
        for job in metadata.job:
            indices = objects.index[(objects.job == job) & (objects.object_type == 'PLA')][:10 if job != first_job else 9]
            objects.loc[indices, 'area'] = .02
        result = cc.comparison_results(objects, metadata, min_per_stratum=10)
        arm = next(iter(ARMS.values()))
        weight = result['standard_weights'].query('pla_arm == @arm and size_index == 0').iloc[0]
        self.assertFalse(weight.included)
        self.assertEqual(weight.min_animal_n, 9)
        self.assertIn(ANIMAL_ORDER[0], weight.reason)
        selected = result['standardization'].query('pla_arm == @arm and estimator == "restricted"').iloc[0]
        self.assertEqual(selected.min_retained_pct, 90.)
        self.assertEqual(selected.max_retained_pct, 91.)
        self.assertIn('<0.03', selected.excluded_strata)

    def test_actual_overlap_area_is_not_area_of_positive_objects(self):
        objects, metadata = fixture([50] * 6)
        result = cc.comparison_results(objects, metadata)
        np.testing.assert_allclose(result['animal_metrics'].overlap_area_ratio_pct, 5.)
        np.testing.assert_allclose(result['animal_metrics'].association_pct, 50.)

    def test_missing_values_remain_unknown_and_strict_payload_serializes_null(self):
        objects, metadata = fixture()
        objects['chat_present'] = objects.chat_present.astype(float)
        objects.loc[0, 'chat_present'] = np.nan
        objects.loc[1, 'area'] = np.nan
        result = cc.comparison_results(objects, metadata)
        arm = next(iter(ARMS.values()))
        animal = result['animal_metrics'].query('animal_id == "21F" and pla_arm == @arm').iloc[0]
        self.assertTrue(np.isnan(animal.association_pct))
        self.assertTrue(np.isnan(animal.mean_pla_area_um2))
        self.assertTrue(np.isnan(animal.large_pla_pct))
        contrast = result['contrasts'].query('pla_arm == @arm and metric == "association_pct"').iloc[0]
        self.assertTrue(np.isnan(contrast.home_mean))
        self.assertEqual(contrast.n_finite_home, 2)
        self.assertTrue(np.isnan(result['exact_summary'].query('pla_arm == @arm').two_sided_p.iloc[0]))
        payload = cc.comparison_payload(result)
        json.dumps(payload, allow_nan=False)
        self.assertTrue(any(row['association_pct'] is None for row in payload['animal_metrics']))

    def test_rejects_wrong_groups_repeated_animals_extra_jobs_and_invalid_flags(self):
        objects, metadata = fixture()
        with self.assertRaisesRegex(ValueError, 'twelve unique'):
            cc.comparison_results(objects, pd.concat([metadata, metadata.iloc[[0]]]))
        bad = metadata.copy()
        bad.loc[0, 'running_group'] = cc.RUNNING
        with self.assertRaisesRegex(ValueError, 'assignments'):
            cc.comparison_results(objects, bad)
        extra = objects.iloc[[0]].copy()
        extra['job'] = 'interrupted_duplicate'
        with self.assertRaisesRegex(ValueError, 'primary sections'):
            cc.comparison_results(pd.concat([objects, extra]), metadata)
        bad_objects = objects.copy()
        bad_objects['chat_present'] = bad_objects.chat_present.astype(object)
        bad_objects.loc[0, 'chat_present'] = 'false'
        with self.assertRaisesRegex(ValueError, 'binary'):
            cc.comparison_results(bad_objects, metadata)
        for threshold in [0, True, 1.5]:
            with self.assertRaisesRegex(ValueError, 'positive integer'):
                cc.comparison_results(objects, metadata, min_per_stratum=threshold)


    def reviewed_fixture(self):
        objects, metadata = fixture()
        objects['x'] = 10.
        objects['y'] = 10.
        records = metadata.to_dict('records')
        for record in records:
            record.update(reviewed_polygon_native=[[-1, -1], [1, -1], [1, 1], [-1, 1]],
                          landmark_reviewed=True, orientation='none', orientation_reviewed=True,
                          scale_common_confirmed=True, support_reviewed=True, support_qc='complete',
                          reviewer='Synthetic test reviewer', second_reviewer='')
        return objects, records

    def test_registration_gate_rejects_empty_partial_unready_and_duplicate_reviews(self):
        objects, records = self.reviewed_fixture()
        with mock.patch('bio_visuals.handoff.registered_fraction_maps') as register:
            for subset in ([], records[:-1]):
                status = cc.condition_registration(objects, subset)
                self.assertFalse(status['ready'])
                self.assertEqual(status['maps'], {})
                self.assertEqual(status['total'], len(subset))
            partial = [dict(record) for record in records]
            partial[0]['support_reviewed'] = False
            status = cc.condition_registration(objects, partial)
            self.assertFalse(status['ready'])
            self.assertEqual(status['approved'], 11)
            self.assertEqual(len(status['needs_review']), 1)
            self.assertIn('support', status['needs_review'][0]['needs_review'])
            duplicate = records[:-1] + records[:1]
            status = cc.condition_registration(objects, duplicate)
            self.assertFalse(status['ready'])
            self.assertIn('identities', status['message'])
            register.assert_not_called()

    def test_registration_adapter_reverses_direction_and_preserves_unknown_bins(self):
        objects, records = self.reviewed_fixture()
        expected_map = {'Home cage': dict(mean=np.array([[60., np.nan]]), n=np.array([[3, 1]])),
                        'Running': dict(mean=np.array([[40., np.nan]]), n=np.array([[3, 2]])),
                        'difference': np.array([[-20., np.nan]])}
        source = dict(x_edges=np.array([0., 10., 20.]), y_edges=np.array([0., 10.]),
                      arms={arm: expected_map for arm in ARMS.values()})
        with mock.patch('bio_visuals.handoff.registered_fraction_maps', return_value=source) as register:
            status = cc.condition_registration(objects, records, width=10., min_pla=5)
            register.assert_called_once_with(objects, records, width=10., min_pla=5)
        self.assertTrue(status['ready'])
        self.assertEqual(status['approved'], 12)
        self.assertEqual(status['needs_review'], [])
        for result in status['maps'].values():
            self.assertEqual(result['home_minus_running'], [[20., None]])
            self.assertEqual(result['home_mean'], [[60., None]])
            self.assertEqual(result['home_n'], [[3, 1]])
            self.assertEqual(result['running_n'], [[3, 2]])
        json.dumps(status, allow_nan=False)

    def test_reviewed_registration_uses_equal_animals_and_observed_common_support(self):
        objects, records = self.reviewed_fixture()
        status = cc.condition_registration(objects, records)
        self.assertTrue(status['ready'])
        for result in status['maps'].values():
            difference = np.asarray(result['home_minus_running'], float)
            np.testing.assert_allclose(difference[np.isfinite(difference)], [80.])
            self.assertEqual(np.isfinite(difference).sum(), 1)
            self.assertEqual(np.asarray(result['home_n']).max(), 3)


if __name__ == '__main__':
    unittest.main()
