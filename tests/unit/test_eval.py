""" Unit tests for the Eval class."""

from pdm_utils.classes import eval
import unittest

class TestEvalClass(unittest.TestCase):


    def setUp(self):
        self.eval_result = eval.Eval()

    def test_default_eval_1(self):
        """Verify default status is empty."""
        self.assertEqual(self.eval_result.status, "")



if __name__ == '__main__':
    unittest.main()
