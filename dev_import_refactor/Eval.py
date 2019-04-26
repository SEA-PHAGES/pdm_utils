"""Represents a structure to contain results of an evaluation.
"""






class EvalResult:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.status = '' #correct, warning, error, etc.
        self.messages = {\
        "correct":"empty",\
        "warning":"empty",\
        "error":"empty"} #Keys should match the possible status options


    def current_message(self):

        try:
            message = self.messages[self.status]
        except:
            message = "No matching message available."

        return message


























###
