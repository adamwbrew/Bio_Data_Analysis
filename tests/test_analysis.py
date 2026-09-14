"""Small scientific invariants independent of the large, untracked input CSVs.

Run from the repository root with ``python -m unittest discover -s tests -v``.
These intentionally unequal synthetic objects expose denominator confusion,
tie handling, and spatial support mistakes that real-data smoke tests can miss.
"""

import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from bio_visuals.data import REGIONS, quality_control, summarize
from bio_visuals.plots import area_sensitivity, ecdf_figures, overview
from bio_visuals.spatial import spatial_grid_data


def synthetic_objects():
    """Four jobs, each with four PLA objects and two ChAT objects.

    PLA areas are 1, 1, 2, 8 and overlap areas 1, 0, 0.5, 0. Therefore
    object association = 50%, mean object overlap = 31.25%, pooled area
    overlap = 12.5%. One ChAT object also has PLA present and must not be
    counted as a fifth PLA object. Only the first spatial bin has two PLA.
    """
    rows = []
    for job_index, (job, region) in enumerate(REGIONS.items()):
        shift = 100.0 * job_index
        specs = [
            ("PLA", 1.0, 1.0, 0.1, True, True, 1.0),
            ("PLA", 1.0, 0.0, 0.1, True, False, 3.0),
            ("PLA", 2.0, 0.5, 0.2, True, True, 11.0),
            ("PLA", 8.0, 0.0, 0.8, True, False, 21.0),
            ("ChAT", 4.0, 2.0, 0.02, True, True, 1.0),
            ("ChAT", 10.0, 0.0, 0.03, False, True, 31.0),
        ]
        for object_id, (typ, area, overlap, intensity, pla, chat, x) in enumerate(specs):
            x += shift
            y = shift + 1.0
            rows.append({
                "job": job, "region": region, "object_id": object_id,
                "object_type": typ, "area": area, "inner_area": 0.0,
                "outer_area": area, "overlap_area": overlap,
                "overlap_pct": 100.0 * overlap / area,
                "intensity": intensity, "pla_present": pla,
                "chat_present": chat,
                "min_diameter": 0.0 if object_id == 0 else 1.0,
                "median_diameter": 1.0, "max_diameter": 2.0,
                "xmin": x if object_id == 0 else x - 0.5,
                "xmax": x if object_id == 0 else x + 0.5,
                "ymin": y - 0.5, "ymax": y + 0.5,
                "x": x, "y": y,
            })
    return pd.DataFrame(rows)


class AnalysisInvariantTests(unittest.TestCase):
    def setUp(self):
        self.objects = synthetic_objects()

    def tearDown(self):
        plt.close("all")

    def test_summary_keeps_object_types_and_denominators_distinct(self):
        table = summarize(self.objects)
        self.assertEqual(len(table), 8)
        self.assertEqual(table.n_objects.sum(), 24)
        self.assertEqual(int(self.objects.pla_present.sum()), 20)
        self.assertEqual(int(self.objects.object_type.eq("PLA").sum()), 16)
        for job in REGIONS:
            pla = table[(table.job == job) & (table.object_type == "PLA")].iloc[0]
            chat = table[(table.job == job) & (table.object_type == "ChAT")].iloc[0]
            self.assertEqual(pla.n_objects, 4)
            self.assertEqual(chat.n_objects, 2)
            self.assertEqual(pla.other_marker_count, 2)
            self.assertAlmostEqual(pla.other_marker_pct, 50.0)
            self.assertAlmostEqual(pla.mean_object_overlap_pct, 31.25)
            self.assertAlmostEqual(pla.area_weighted_overlap_pct, 12.5)
            self.assertAlmostEqual(chat.area_weighted_overlap_pct, 100.0 * 2.0 / 14.0)
            self.assertAlmostEqual(pla.median_area_um2, 1.5)

    def test_overview_uses_pla_object_type_in_all_three_measures(self):
        _, table = overview(self.objects)
        np.testing.assert_array_equal(table.pla_objects, [4] * 4)
        np.testing.assert_array_equal(table.chat_objects, [2] * 4)
        np.testing.assert_array_equal(table.pla_with_chat, [2] * 4)
        np.testing.assert_allclose(table.pla_object_fraction_pct, 50.0)
        np.testing.assert_allclose(table.pla_mean_overlap_pct, 31.25)
        np.testing.assert_allclose(table.pla_area_weighted_overlap_pct, 12.5)

    def test_ecdf_accumulates_ties_and_reaches_one_hundred(self):
        before = self.objects.copy(deep=True)
        _, table = ecdf_figures(self.objects)
        selected = table[(table.job == "4594") & (table.object_type == "PLA") & (table.metric == "area")]
        np.testing.assert_array_equal(selected.value, [1.0, 2.0, 8.0])
        np.testing.assert_array_equal(selected.object_count_at_value, [2, 1, 1])
        np.testing.assert_allclose(selected.cumulative_pct, [50.0, 75.0, 100.0])
        for (_, typ, _), group in table.groupby(["job", "object_type", "metric"]):
            self.assertEqual(group.object_count_at_value.sum(), 4 if typ == "PLA" else 2)
            self.assertAlmostEqual(group.cumulative_pct.iloc[-1], 100.0)
            self.assertTrue(group.cumulative_pct.is_monotonic_increasing)
        assert_frame_equal(self.objects, before)

    def test_area_sensitivity_is_inclusive_and_preserves_raw_data(self):
        before = self.objects.copy(deep=True)
        _, table = area_sensitivity(self.objects, thresholds=[0.0, 1.0, 2.0, 8.0, 9.0])
        for _, group in table.groupby("job"):
            np.testing.assert_array_equal(group.retained_pla_count, [4, 4, 2, 1, 0])
            np.testing.assert_allclose(group.retained_pla_pct, [100, 100, 50, 25, 0])
            np.testing.assert_array_equal(group.pla_with_chat_count, [2, 2, 1, 0, 0])
            np.testing.assert_allclose(group.pla_with_chat_pct.iloc[:4], [50, 50, 50, 0])
            self.assertTrue(np.isnan(group.pla_with_chat_pct.iloc[-1]))
        assert_frame_equal(self.objects, before)

    def test_qc_retains_zero_geometry_and_detects_deliberate_inconsistency(self):
        before = self.objects.copy(deep=True)
        qc = quality_control(self.objects)
        np.testing.assert_array_equal(qc.any_zero_diameter, [1] * 4)
        np.testing.assert_array_equal(qc.zero_bbox_width_or_height, [1] * 4)
        np.testing.assert_array_equal(qc.overlap_presence_mismatch, [0] * 4)
        np.testing.assert_array_equal(qc.overlap_pct_formula_mismatch, [0] * 4)
        # Make one export flag and one exported percentage disagree with area.
        inconsistent = self.objects.copy(deep=True)
        inconsistent.loc[0, "chat_present"] = False
        inconsistent.loc[0, "overlap_pct"] = 75.0
        flagged = quality_control(inconsistent).set_index("job")
        self.assertEqual(flagged.loc["4594", "overlap_presence_mismatch"], 1)
        self.assertEqual(flagged.loc["4594", "overlap_pct_formula_mismatch"], 1)
        self.assertEqual(flagged.loc["4594", "n_objects"], 6)
        assert_frame_equal(self.objects, before)

    def test_spatial_bins_conserve_objects_and_use_raw_local_support(self):
        before = self.objects.copy(deep=True)
        grid = spatial_grid_data(self.objects, bin_width=10, min_support=2, metric="pla_fraction")
        for _, group in grid.groupby("job"):
            self.assertEqual(group.pla_count.sum(), 4)
            self.assertEqual(group.chat_count.sum(), 2)
            self.assertEqual(group.pla_with_chat_count.sum(), 2)
            self.assertEqual(group.raw_pla_support_pass.sum(), 1)
            supported = group[group.raw_pla_support_pass].iloc[0]
            self.assertEqual(supported.pla_count, 2)
            self.assertAlmostEqual(supported.display_value, 50.0)
            self.assertTrue(group.loc[~group.raw_pla_support_pass, "display_value"].isna().all())
            self.assertFalse(group.roi_membership_known.any())
        assert_frame_equal(self.objects, before)

    def test_spatial_smoothing_uses_ratio_of_weighted_counts_without_changing_raw_counts(self):
        kwargs = dict(bin_width=10, min_support=2, metric="pla_fraction")
        raw = spatial_grid_data(self.objects, smoothing=0, **kwargs)
        smooth = spatial_grid_data(self.objects, smoothing=1, **kwargs)
        raw_columns = ["job", "x_left", "y_top", "pla_count", "chat_count", "pla_with_chat_count", "pla_fraction_pct", "mean_overlap_pct", "raw_pla_support_pass"]
        assert_frame_equal(raw[raw_columns], smooth[raw_columns])
        # First bin has [2, 1, 1, 0] PLA and [1, 1, 0, 0] associations
        # at distances 0,1,2,3. Gaussian edge normalization cancels in ratio.
        expected = 100 * (1 + np.exp(-0.5)) / (2 + np.exp(-0.5) + np.exp(-2.0))
        np.testing.assert_allclose(smooth.loc[smooth.raw_pla_support_pass, "display_value"], expected)
        self.assertTrue(smooth.loc[~smooth.raw_pla_support_pass, "display_value"].isna().all())
        count_display = spatial_grid_data(self.objects, bin_width=10, min_support=2, smoothing=1, metric="pla_count")
        self.assertTrue(count_display.loc[count_display.pla_count == 0, "display_value"].isna().all())

    def test_spatial_invalid_parameters_fail_clearly(self):
        for kwargs in [
            {"bin_width": 0}, {"bin_width": -10}, {"min_support": 0},
            {"min_support": 1.5}, {"smoothing": -1}, {"metric": "tissue_density"},
        ]:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    spatial_grid_data(self.objects, **kwargs)


if __name__ == "__main__":
    unittest.main()
