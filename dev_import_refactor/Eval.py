"""Represents a structure to contain results of an evaluation.
"""






class EvalResult:

    # Initialize all attributes:
    def __init__(self):

        # Initialize all non-calculated attributes:
        self.status = "" #correct, warning, error, etc.
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















def construct_warning(message_warning,message_error):

    eval = EvalResult()
    eval.status = "warning"
    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval


def construct_error(message_error):

    eval = EvalResult()
    eval.status = "error"
    eval.messages["error"] = message_error

    return eval







###
