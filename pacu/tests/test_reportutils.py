import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pacu.app.utils import reportutils


class TestReportUtils(unittest.TestCase):
    """
    Tests for the report utils.
    """

    def test_create_upsetplot_overlap(self) -> None:
        """
        Tests the 'create_upsetplot_overlap' function.
        :return: None
        """
        data_overlap = pd.DataFrame([
            {'chr': 'c1', 'pos': 1, 'phages': True, 'gubbins': False, 'depth': False},
            {'chr': 'c1', 'pos': 2, 'phages': True, 'gubbins': False, 'depth': True},
        ])
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_png = Path(dir_, 'overlap.png')
            reportutils.create_upsetplot_overlap(data_overlap, path_png)
            self.assertGreater(path_png.stat().st_size, 0)

    def test_create_upsetplot_overlap_empty(self) -> None:
        """
        Tests the 'create_upsetplot_overlap' function when no positions are filtered.
        :return: None
        """
        data_overlap = pd.DataFrame(columns=['chr', 'pos', 'phages', 'gubbins', 'depth'])
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_png = Path(dir_, 'overlap.png')
            reportutils.create_upsetplot_overlap(data_overlap, path_png)
            self.assertGreater(path_png.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
