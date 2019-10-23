"""Unit tests for main run script."""

import unittest
from unittest.mock import patch
import argparse
from pdm_utils import run

class TestRunFunctions1(unittest.TestCase):

    def setUp(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("pipeline")
        self.args1 = self.parser.parse_args(["import_dev"])

    @patch("pdm_utils.pipelines.db_import.import_genome.main")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_1(self, parse_args_mock, run_pipeline_mock):
        """Verify that main runs correctly with:
        valid pipeline and valid default options."""
        parse_args_mock.return_value = self.args1
        run.main()
        with self.subTest():
            self.assertTrue(parse_args_mock.called)
        with self.subTest():
            self.assertTrue(run_pipeline_mock.called)



if __name__ == '__main__':
    unittest.main()
