"""Protect missing estimates, support masks, and offline embedding in ML views."""
import json
import re
import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bio_visuals import morphology_views as views


class MorphologyViewsTests(unittest.TestCase):
    def test_missing_estimates_are_null_not_zero_or_nonstandard_json(self):
        result = {'animals': pd.DataFrame([{'pla_arm': 'TDP43–HK1',
                                           'association_pct': np.nan,
                                           'proportion_pct': 0.}]),
                  'profiles': pd.DataFrame([{'zero_min_pct': np.inf}])}
        encoded = json.dumps(views.explorer_payload(result), allow_nan=False)
        decoded = json.loads(encoded)
        self.assertIsNone(decoded['animals'][0]['association_pct'])
        self.assertEqual(decoded['animals'][0]['proportion_pct'], 0.)
        self.assertIsNone(decoded['profiles'][0]['zero_min_pct'])

    def test_native_fraction_masks_sparse_bins_but_keeps_supported_zero(self):
        section = dict(job='a', animal_id='21F', running_group='Home cage',
                       pla_arm='TDP43–HK1', extent=[0, 1000, 0, 1000], min_count=10,
                       n=[[9, 10], [0, 12]], n_all=[[10, 10], [5, 15]],
                       cluster_fractions={'1': [[40., 0.], [100., 25.]]})
        fig = views.spatial_figure({'spatial': [section]}, 'TDP43–HK1')
        try:
            shown = fig.axes[0].images[1].get_array()
            self.assertTrue(shown.mask[0, 0])
            self.assertTrue(shown.mask[1, 0])
            self.assertFalse(shown.mask[0, 1])
            self.assertEqual(shown[0, 1], 0.)
            self.assertEqual(shown[1, 1], 25.)
            # An excluded-only bin is observed PLA support, not absent tissue.
            support = fig.axes[0].images[0].get_array()
            self.assertFalse(np.ma.getmaskarray(support)[1, 0])
            self.assertEqual(section['cluster_fractions']['1'][0][0], 40.)
        finally:
            plt.close(fig)

    def test_fragment_ids_do_not_collide_and_script_text_is_escaped(self):
        result = {'definitions': {'source': '</script><script>alert(1)</script>'}}
        first = views.explorer_html(result, include_plotly=False)
        second = views.explorer_html(result, include_plotly=False)
        first_id = re.search(r'id="(morphology-[a-z0-9]+)"', first).group(1)
        second_id = re.search(r'id="(morphology-[a-z0-9]+)"', second).group(1)
        self.assertNotEqual(first_id, second_id)
        self.assertNotIn('</script><script>alert(1)</script>', first)
        self.assertIn('\\u003c/script>', first)
        self.assertNotIn('__ROOT_ID__', first)
        self.assertNotIn('__DATA__', first)


if __name__ == '__main__':
    unittest.main()
