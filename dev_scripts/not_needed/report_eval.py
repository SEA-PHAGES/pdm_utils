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




def eval_status(phage_name,import_status):

    eval = Eval.EvalResult()
    eval.status = "warning"
    message_warning = \
    "The status %s is not currently in the database." % import_status

    message_error = \
    "The status %s is not correct for %s." % (import_status,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval




def eval_description_field(phage_name,description_field):

    eval = Eval.EvalResult()
    eval.status = "warning"
    message_warning = \
    "The description field %s is not commonly used." % description_field

    message_error = \
    "The description field %s is not correct for %s."
    % (description_field,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval






def eval_field_not_valid(phage_name,field):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = \
    "The %s field is not valid for %s." \
    % (field,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval










def eval_phageid_not_present(phage_name):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = \
    "The phage %s is not in the database." \
    % phage_name

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval




def eval_phageid_present(phage_name):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = \
    "The phage %s is already in the database." \
    % phage_name

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval







def eval_field_empty(phage_name,field):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = \
    "The %s field for phage %s should not be empty."
    #"The phage %s does not have correctly populated %s field." \
    % (field,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval





def eval_field_not_empty(phage_name,field):

    eval = Eval.EvalResult()
    eval.status = "error"
    message_warning = "none"
    message_error = \
    "The %s field for phage %s should be empty."
    % (field,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval





def eval_status_discrepancy(phage_name,field):

    eval = Eval.EvalResult()
    eval.status = "warning"
    message_warning = \
    "The phage %s to be added is listed as Final status, but no Draft \
    (or other) genome is listed to be removed."
    % field

    message_error = \
    "The status is incorrect for %s."
    % (field,phage_name)

    eval.messages["warning"] = message_warning
    eval.messages["error"] = message_error

    return eval
















###CDS feature









###tRNA feature






###Genome










###Ticket compared to Test genome






###Test genome compared to PhagesDB genome












###Test genome compare to Phamerator genome








###
