"""Represents a structure to contain results of an evaluation.
"""






class Eval:

    def __init__(self, id = "", definition = "", status = "", result = ""):

        self.id = id # Unique identifier for the specific evaluation.
        self.definition = definition # Description of what was evaluated.
        self.status = status # Not_evaluated, Ignored, Correct, Warning, Error.
        self.result = result # Customized message reporting details of
                         # the evaluation.






#
# def construct_warning(message_warning, message):
#
#     eval = EvalResult()
#     eval.status = "warning"
#     eval.result = message
#
#     return eval
#
#
# def construct_error(message_error):
#
#     eval = EvalResult()
#     eval.status = "error"
#     eval.messages["error"] = message_error
#
#     return eval
#
#
#
# def construct_other(status_other, message_other):
#
#     eval = EvalResult()
#     eval.status = status_other
#     eval.messages[status_other] = message_other
#
#     return eval
#










###
