""" Unit tests for the Eval class."""



import Eval
import unittest




class TestEvalClass(unittest.TestCase):


    def setUp(self):
        self.eval_result = Eval.EvalResult()



    def test_default_eval_1(self):
        """Verify default status is empty."""
        self.assertEqual(self.eval_result.status, "")

    def test_default_eval_2(self):
        """Verify default message dictionary length."""
        self.assertEqual(len(self.eval_result.messages), 3)



    def test_current_message_1(self):
        """Verify default retrieved message."""
        message = "No matching message available."
        self.assertEqual(self.eval_result.current_message(), message)

    def test_current_message_2(self):
        """Verify retrieved correct message."""
        self.eval_result.status = "correct"
        self.eval_result.messages["correct"] = "empty_correct"
        self.assertEqual(self.eval_result.current_message(), "empty_correct")

    def test_current_message_3(self):
        """Verify retrieved warning message."""
        self.eval_result.status = "warning"
        self.eval_result.messages["warning"] = "empty_warning"
        self.assertEqual(self.eval_result.current_message(), "empty_warning")

    def test_current_message_4(self):
        """Verify retrieved error message."""
        self.eval_result.status = "error"
        self.eval_result.messages["error"] = "empty_error"
        self.assertEqual(self.eval_result.current_message(), "empty_error")

    def test_current_message_5(self):
        """Verify retrieved non-standard message."""
        self.eval_result.status = "other"
        self.eval_result.messages["other"] = "empty_other"
        self.assertEqual(self.eval_result.current_message(), "empty_other")




    def test_construct_warning_1(self):
        """Verify construction of warning result."""
        self.warning_result = Eval.construct_warning("warning_message", "error_message")
        with self.subTest():
            self.assertEqual(self.warning_result.status, "warning")
        with self.subTest():
            self.assertEqual(self.warning_result.messages["warning"], "warning_message")
        with self.subTest():
            self.assertEqual(self.warning_result.messages["error"], "error_message")
        with self.subTest():
            self.assertEqual(self.warning_result.current_message(), "warning_message")


    def test_construct_error_1(self):
        """Verify construction of error result."""
        self.error_result = Eval.construct_error("error_message")
        with self.subTest():
            self.assertEqual(self.error_result.status, "error")
        with self.subTest():
            self.assertEqual(self.error_result.messages["warning"], "empty")
        with self.subTest():
            self.assertEqual(self.error_result.messages["error"], "error_message")
        with self.subTest():
            self.assertEqual(self.error_result.current_message(), "error_message")


    def test_construct_other_1(self):
        """Verify construction of other result."""
        self.other_result = Eval.construct_other("other","other_message")
        with self.subTest():
            self.assertEqual(self.other_result.status, "other")
        with self.subTest():
            self.assertEqual(self.other_result.messages["warning"], "empty")
        with self.subTest():
            self.assertEqual(self.other_result.messages["error"], "empty")
        with self.subTest():
            self.assertEqual(self.other_result.messages["other"], "other_message")
        with self.subTest():
            self.assertEqual(self.other_result.current_message(), "other_message")









if __name__ == '__main__':
    unittest.main()
