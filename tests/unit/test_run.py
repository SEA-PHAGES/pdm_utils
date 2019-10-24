"""Unit tests for main run script."""

import unittest
from unittest.mock import patch
from pdm_utils import run

class TestRunFunctions1(unittest.TestCase):

    @patch("pdm_utils.pipelines.db_import.import_genome.main")
    def test_main_1(self, run_pipeline_mock):
        """Verify that main runs correctly with:
        valid pipeline and valid default options."""
        unparsed_args = ["pdm_utils.run", "import_dev"]
        run.main(unparsed_args)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)



if __name__ == '__main__':
    unittest.main()
