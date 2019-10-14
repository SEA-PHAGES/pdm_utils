"""Unit tests for constants."""

import unittest
from pdm_utils.pipelines.db_import import run_mode
from unittest.mock import patch

class TestConstants(unittest.TestCase):

    def test_get_eval_flag_dict_1(self):
        """Verify 'base' eval dictionary is returned."""
        dict = run_mode.get_eval_flag_dict("base")
        with self.subTest():
            self.assertEqual(len(dict.keys()), 11)
        with self.subTest():
            self.assertTrue(dict["check_locus_tag"])

    def test_get_eval_flag_dict_2(self):
        """Verify 'pecaan' eval dictionary is returned."""
        dict = run_mode.get_eval_flag_dict("pecaan")
        self.assertFalse(dict["check_locus_tag"])

    def test_get_eval_flag_dict_3(self):
        """Verify 'phagesdb' eval dictionary is returned."""
        dict = run_mode.get_eval_flag_dict("phagesdb")
        self.assertFalse(dict["import_locus_tag"])

    def test_get_eval_flag_dict_4(self):
        """Verify 'sea_auto' eval dictionary is returned."""
        dict = run_mode.get_eval_flag_dict("sea_auto")
        self.assertFalse(dict["check_locus_tag"])

    def test_get_eval_flag_dict_5(self):
        """Verify 'misc' eval dictionary is returned."""
        dict = run_mode.get_eval_flag_dict("misc")
        self.assertFalse(dict["check_locus_tag"])


    @patch("pdm_utils.functions.basic.ask_yes_no", return_value = False)
    def test_get_eval_flag_dict_6(self, mock_ask_yes_no):
        """Verify 'custom' eval dictionary is returned."""
        dict = run_mode.get_eval_flag_dict("custom")
        with self.subTest():
            self.assertEqual(len(dict.keys()), 11)
        with self.subTest():
            self.assertFalse(dict["check_locus_tag"])


###
