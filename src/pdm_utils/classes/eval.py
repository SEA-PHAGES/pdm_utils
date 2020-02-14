"""Represents a structure to contain results of an evaluation.
"""

class Eval:
    def __init__(self, id="", definition="", result="", status=""):
        self.id = id # Unique identifier for the specific evaluation.
        self.definition = definition # Description of what was evaluated.
        self.status = status # untested, ignored, correct, warning, error.
        self.result = result # Customized message reporting details of
                             # the evaluation.

    def __str__(self):
        return (f"Evaluation ID: {self.id}. "
                f"Status: {self.status}. "
                f"Definition: {self.definition} "
                f"Result: {self.result}")
