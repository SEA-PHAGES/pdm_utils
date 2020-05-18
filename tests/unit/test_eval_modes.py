"""Unit tests for eval_modes."""

import unittest
from unittest.mock import patch

from pdm_utils.functions import eval_modes

class TestRunModes(unittest.TestCase):

    def test_get_eval_flag_dict_1(self):
        """Verify 'base' evaluation dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("base")
        with self.subTest():
            self.assertEqual(len(dict.keys()), 14)
        with self.subTest():
            self.assertTrue(dict["check_locus_tag"])

    def test_get_eval_flag_dict_2(self):
        """Verify 'draft' evaluation dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("draft")
        self.assertFalse(dict["check_locus_tag"])

    def test_get_eval_flag_dict_3(self):
        """Verify 'final' evaluation dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("final")
        self.assertFalse(dict["import_locus_tag"])

    def test_get_eval_flag_dict_4(self):
        """Verify 'auto' evaluation dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("auto")
        self.assertFalse(dict["check_locus_tag"])

    def test_get_eval_flag_dict_5(self):
        """Verify 'misc' evaluation dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("misc")
        self.assertFalse(dict["check_locus_tag"])


    @patch("pdm_utils.functions.basic.ask_yes_no", return_value = False)
    def test_get_eval_flag_dict_6(self, mock_ask_yes_no):
        """Verify 'custom' evaluation dictionary is returned."""
        dict = eval_modes.get_eval_flag_dict("custom")
        with self.subTest():
            self.assertEqual(len(dict.keys()), 14)
        with self.subTest():
            self.assertFalse(dict["check_locus_tag"])

if __name__ == '__main__':
    unittest.main()
