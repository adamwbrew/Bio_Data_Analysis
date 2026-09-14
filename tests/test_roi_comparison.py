"""Scientific and loading invariants for repeat ROI analyses."""

import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from bio_visuals.data import (load_data, ALL_REGIONS, REGIONS, ROI_PAIRS,
                             SOURCE_COLUMNS, FILE_PREFIX, FILE_SUFFIX)
from bio_visuals.comparison import (comparison_summary, exact_signature_comparison,
                                    paired_grid_data, paired_spatial_figure)
from bio_visuals.plots import overview
from bio_visuals.spatial import coordinate_overview
from test_analysis import synthetic_objects


def paired_objects():
    box = synthetic_objects()
    box["roi_mode"] = "box"
    box["pair_id"] = box.region.str.replace(" PLA", "", regex=False)
    outline = box.copy(deep=True)
    outline["job"] = outline.job.map(ROI_PAIRS)
    outline["region"] = outline.job.map(ALL_REGIONS)
    outline["roi_mode"] = "outline"
    # Remove the large, nonassociated PLA object in each repeat.
    outline = outline.loc[outline.object_id != 3].copy()
    # A different ID for the same measurements must still match.
    outline["object_id"] += 1000
    return pd.concat([box, outline], ignore_index=True)


class ROIComparisonTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_paired_metrics_keep_runs_separate_and_report_percentage_points(self):
        objects = paired_objects()
        before = objects.copy(deep=True)
        table = comparison_summary(objects)
        pla = table.loc[table.object_type == "PLA"]
        np.testing.assert_array_equal(pla.box_count, [4] * 4)
        np.testing.assert_array_equal(pla.outline_count, [3] * 4)
        np.testing.assert_allclose(pla.count_delta_pct, -25)
        np.testing.assert_allclose(pla.association_delta, 100 * 2 / 3 - 50)
        np.testing.assert_allclose(pla.weighted_overlap_delta, 37.5 - 12.5)
        assert_frame_equal(objects, before)

    def test_signatures_ignore_ids_and_count_duplicates_without_many_to_many_join(self):
        objects = paired_objects()
        dup = objects.loc[(objects.job == "4594") & (objects.object_id == 0)].copy()
        dup["object_id"] = 999
        objects = pd.concat([objects, dup], ignore_index=True)
        table = exact_signature_comparison(objects)
        row = table.loc[(table.pair_id == "22M_Pfkp-TDP43") & (table.object_type == "PLA")].iloc[0]
        self.assertEqual(row.identical_measurement_rows, 3)
        self.assertEqual(row.box_rows_without_identical_signature, 2)
        self.assertEqual(row.outline_rows_without_identical_signature, 0)
        # Changing a measurement with unchanged geometry must not count as exact.
        objects.loc[(objects.job == "4999") & (objects.object_id == 1000), "intensity"] += 0.01
        changed = exact_signature_comparison(objects)
        row = changed.loc[(changed.pair_id == "22M_Pfkp-TDP43") & (changed.object_type == "PLA")].iloc[0]
        self.assertEqual(row.identical_measurement_rows, 2)
        self.assertEqual(row.outline_rows_without_identical_signature, 1)

    def test_aligned_grid_differences_conserve_counts_in_both_runs(self):
        grid = paired_grid_data(paired_objects(), bin_width=10)
        for (pair, typ), group in grid.groupby(["pair_id", "object_type"]):
            self.assertEqual(group.box_count.sum(), 4 if typ == "PLA" else 2)
            self.assertEqual(group.outline_count.sum(), 3 if typ == "PLA" else 2)
            self.assertEqual(group.count_delta.sum(), -1 if typ == "PLA" else 0)
            np.testing.assert_array_equal(group.count_delta, group.outline_count - group.box_count)
            self.assertFalse(group.duplicated(["x_left", "y_top"]).any())
        fig = paired_spatial_figure(grid, "22M_Pfkp-TDP43")
        self.assertEqual(len(fig.axes), 12)  # six panels and six colorbars

    def test_main_views_reject_pooled_repeat_runs(self):
        objects = paired_objects()
        for plot in (overview, coordinate_overview):
            with self.assertRaisesRegex(ValueError, "Select one roi_mode"):
                plot(objects)

    def test_loader_selects_known_roi_modes_without_requiring_unselected_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for job, region in ALL_REGIONS.items():
                row = ['H:\\source.ome.tiff', region, 'algorithm', 0, 'PLA - TRITC_E1',
                       1, 100, 1, 0, 1, 0.2, 0, 1, 0.5, 1, 1, 1, 2, 3, 4]
                with (root / (FILE_PREFIX + job + FILE_SUFFIX)).open('w', encoding='utf-8-sig', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(SOURCE_COLUMNS)
                    writer.writerow(row)
            with patch('bio_visuals.data.KNOWN_ROWS', {job: 1 for job in ALL_REGIONS}):
                both, manifest = load_data(root, roi_mode='all')
                self.assertEqual(len(both), 8)
                self.assertEqual(manifest.sha256.str.len().unique().tolist(), [64])
                self.assertEqual(both.groupby('pair_id').size().tolist(), [2, 2, 2, 2])
                box, _ = load_data(root)
                outline, _ = load_data(root, roi_mode='outline')
                self.assertEqual(set(box.job), set(REGIONS))
                self.assertEqual(set(outline.job), set(ROI_PAIRS.values()))
                (root / (FILE_PREFIX + '4999' + FILE_SUFFIX)).unlink()
                self.assertEqual(len(load_data(root, roi_mode='box')[0]), 4)
                with self.assertRaisesRegex(FileNotFoundError, 'Missing'):
                    load_data(root, roi_mode='outline')
            with self.assertRaisesRegex(ValueError, 'roi_mode'):
                load_data(root, roi_mode='pooled')


if __name__ == '__main__':
    unittest.main()
