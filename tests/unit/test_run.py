"""Unit tests for main run script."""

import unittest
from unittest.mock import patch
from pdm_utils import run

class TestRunFunctions1(unittest.TestCase):

    @patch("pdm_utils.pipelines.get_data.main")
    def test_main_1(self, pipeline_mock):
        """Verify that get_data pipeline is called."""
        unparsed_args = ["pdm_utils.run", "get_data"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.get_db.main")
    def test_main_2(self, pipeline_mock):
        """Verify that get_db pipeline is called."""
        unparsed_args = ["pdm_utils.run", "get_db"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.update_field.main")
    def test_main_3(self, pipeline_mock):
        """Verify that update pipeline is called."""
        unparsed_args = ["pdm_utils.run", "update"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.find_domains.main")
    def test_main_4(self, pipeline_mock):
        """Verify that find_domains pipeline is called."""
        unparsed_args = ["pdm_utils.run", "find_domains"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.import_genome.main")
    def test_main_5(self, pipeline_mock):
        """Verify that import pipeline is called."""
        unparsed_args = ["pdm_utils.run", "import"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.phamerate.main")
    def test_main_6(self, pipeline_mock):
        """Verify that phamerate pipeline is called."""
        unparsed_args = ["pdm_utils.run", "phamerate"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.export_db.main")
    def test_main_7(self, pipeline_mock):
        """Verify that export pipeline is called."""
        unparsed_args = ["pdm_utils.run", "export"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.push_db.main")
    def test_main_8(self, pipeline_mock):
        """Verify that push pipeline is called."""
        unparsed_args = ["pdm_utils.run", "push"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.freeze_db.main")
    def test_main_9(self, pipeline_mock):
        """Verify that freeze pipeline is called."""
        unparsed_args = ["pdm_utils.run", "freeze"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.compare_db.main")
    def test_main_10(self, pipeline_mock):
        """Verify that compare pipeline is called."""
        unparsed_args = ["pdm_utils.run", "compare"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.get_gb_records.main")
    def test_main_11(self, pipeline_mock):
        """Verify that get_gb_records pipeline is called."""
        unparsed_args = ["pdm_utils.run", "get_gb_records"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

    @patch("pdm_utils.pipelines.convert_db.main")
    def test_main_12(self, pipeline_mock):
        """Verify that convert pipeline is called."""
        unparsed_args = ["pdm_utils.run", "convert"]
        run.main(unparsed_args)
        pipeline_mock.assert_called()

if __name__ == '__main__':
    unittest.main()
