""" Unit tests for the Evaluation class."""

import unittest

from pdm_utils.classes import evaluation

class TestEvalClass(unittest.TestCase):


    def setUp(self):
        self.eval_result = evaluation.Evaluation()

    def test_default_eval_1(self):
        """Verify default status is empty."""
        self.assertEqual(self.eval_result.status, "")



if __name__ == '__main__':
    unittest.main()
