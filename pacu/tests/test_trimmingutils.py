import os
import unittest
from pathlib import Path
from unittest import mock

from pacu.app.utils import trimmingutils


class TestTrimmingUtils(unittest.TestCase):
    """
    Tests for the trimming utils.
    """

    def test_trimmomatic_dir_adapters_env_var(self) -> None:
        """
        Tests that the TRIMMOMATIC_ADAPTER_DIR environment variable is used when set.
        :return: None
        """
        with mock.patch.dict(os.environ, {'TRIMMOMATIC_ADAPTER_DIR': '/path/to/adapters', 'CONDA_PREFIX': '/conda'}):
            self.assertEqual(trimmingutils.trimmomatic_dir_adapters(), Path('/path/to/adapters'))

    def test_trimmomatic_dir_adapters_conda(self) -> None:
        """
        Tests that the adapter directory is retrieved from the CONDA prefix.
        :return: None
        """
        with mock.patch.dict(os.environ, {'CONDA_PREFIX': '/conda'}):
            os.environ.pop('TRIMMOMATIC_ADAPTER_DIR', None)
            self.assertEqual(
                trimmingutils.trimmomatic_dir_adapters(), Path('/conda', 'share', 'trimmomatic', 'adapters'))

    def test_trimmomatic_dir_adapters_not_set(self) -> None:
        """
        Tests that an error is raised when the adapter directory cannot be determined.
        :return: None
        """
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(RuntimeError):
                trimmingutils.trimmomatic_dir_adapters()


if __name__ == '__main__':
    unittest.main()
