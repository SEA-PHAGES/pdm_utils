"""Unit tests for eval_modes."""

import unittest
from pdm_utils.functions import eval_modes
from unittest.mock import patch

class TestRunModes(unittest.TestCase):

    def test_get_eval_flag_dict_1(self):
        """Verify 'base' eval dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("base")
        with self.subTest():
            self.assertEqual(len(dict.keys()), 13)
        with self.subTest():
            self.assertTrue(dict["check_locus_tag"])

    def test_get_eval_flag_dict_2(self):
        """Verify 'draft' eval dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("draft")
        self.assertFalse(dict["check_locus_tag"])

    def test_get_eval_flag_dict_3(self):
        """Verify 'final' eval dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("final")
        self.assertFalse(dict["import_locus_tag"])

    def test_get_eval_flag_dict_4(self):
        """Verify 'auto' eval dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("auto")
        self.assertFalse(dict["check_locus_tag"])

    def test_get_eval_flag_dict_5(self):
        """Verify 'misc' eval dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("misc")
        self.assertFalse(dict["check_locus_tag"])


    @patch("pdm_utils.functions.basic.ask_yes_no", return_value = False)
    def test_get_eval_flag_dict_6(self, mock_ask_yes_no):
        """Verify 'custom' eval dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("custom")
        with self.subTest():
            self.assertEqual(len(dict.keys()), 13)
        with self.subTest():
            self.assertFalse(dict["check_locus_tag"])

if __name__ == '__main__':
    unittest.main()
