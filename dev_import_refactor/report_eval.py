import Eval




#generic function



phage_name
def eval_(,,):

    eval = Eval.EvalResult()
    eval.status = "warning"
    message_warning = \
    " %s." % a + \
    "\n %s." % b + \
    "\n %s." % c + \
    "\n ."

    message_error = \
    " %s" % a

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval











###Tickets


def eval_ticket_type(ticket_type):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = "This is not a permissible ticket type: %s" % ticket_type

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval








def eval_field_retrieval(phage_name,field):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = \
    "Unable to retrieve %s data from PhagesDB for %s." \
    % (field,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval




def eval_new_field_value(phage_name,field,new_value):

    eval = Eval.EvalResult()
    eval.status = "warning"
    message_warning = \
    "The following %s is not currently in the database: %s."
    % (field,new_value)

    message_error = \
    "The %s is not correct for: %s."
    % (field,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval









def eval_field_length(field,field_value,length):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = \
    "The following %s designation exceeds the %s character limit: %s." \
    % (field,length,field_value)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval













def eval_data_discrepancy(phage_name,field1,field2):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = \
    "There is a discrepancy between %s and %s data for: %s." \
    % (field1,field2,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval






				write_out(output_file,"\nError: phage %s Cluster and Subcluster discrepancy." % row[1])





###CDS feature









###tRNA feature






###Genome










###Ticket compared to Test genome






###Test genome compared to PhagesDB genome












###Test genome compare to Phamerator genome










#Comparison of host data
def eval_host(phage_name,phamerator_host,import_host):

    eval = Eval.EvalResult()
    eval.status = "warning"
    message_warning = \
    "There is conflicting host data for %s." % phage_name + \
    "\nPhamerator host: %s." % phamerator_host + \
    "\nImport ticket host: %s." % import_host + \
    "\nThe new host data will be imported."

    message_error = \
    "Host data is incorrect for %s" % phage_name

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval



























###
