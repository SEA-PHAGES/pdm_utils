import Eval

#Comparison of host data
def eval_host_check(phage_name,phamerator_host,import_host):

    print("run check")
    eval = Eval.Eval()
    eval.status = "warning"
    message_warning = \
    "There is conflicting host data for %s.\n\
    Phamerator host: %s.\n\
    Import ticket host: %s.\n\
    The new host data will be imported." \
    % (phage_name,phamerator_host,import_host)




    message_error = \
    "Host data is incorrect for %s" % phage_name

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval


def test_print():
    print("hello")






print("test print2")
