"""Represents a structure to contain results of an evaluation.
"""






class Eval:

    def __init__(self, id = "", definition = "", result = "", status = ""):

        self.id = id # Unique identifier for the specific evaluation.
        self.definition = definition # Description of what was evaluated.
        self.status = status # untested, ignored, correct, warning, error.
        self.result = result # Customized message reporting details of
                         # the evaluation.






#
# def construct_warning(message_warning, message):
#
#     evl = EvalResult()
#     evl.status = "warning"
#     evl.result = message
#
#     return evl
#
#
# def construct_error(message_error):
#
#     evl = EvalResult()
#     evl.status = "error"
#     evl.messages["error"] = message_error
#
#     return evl
#
#
#
# def construct_other(status_other, message_other):
#
#     evl = EvalResult()
#     evl.status = status_other
#     evl.messages[status_other] = message_other
#
#     return evl
#










###
