"""Misc. functions to manipulate tickets."""

from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import eval_modes
from pdm_utils.classes import ticket
from pdm_utils.classes import genome

def modify_import_data(data_dict, required_keys, optional_keys, keywords):
    """Modifies ticket data to conform to requirements for an ImportTicket object.

    :param data_dict: Dictionary of import ticket data.
    :type data_dict: dict
    :param required_keys: Set of keys required to be in the data dictionary.
    :type required_keys: set
    :param optional_keys:
        Set of optional keys that are not required to be in the data dictionary.
    :type optional_keys: set
    :param keywords:
        Set of valid keyword values that are handled differently
        than other values.
    :type keywords: set
    :returns: Indicates if the ticket is structured properly.
    :rtype: bool
    """
    # Some import table fields are required. Others are optional.
    # First check if all required fields are present.
    missing_req_fields = required_keys - data_dict.keys()
    # Second, check if there are any fields that are not expected.
    extra_fields = data_dict.keys() - required_keys - optional_keys

    # If row in import table has more fields than headers, they are
    # assigned None by csv.DictReader.
    # If row in import table contains fields with no data, they are
    # assigned "" by csv.DictReader.
    # "retrieve" = get data from PhagesDB.
    # "retain" = get data from current genome in the MySQL database.
    # "parse" = get data from the flatfile itself (only for host_genus and accession)
    # "none" = there is no applicable data (only for accession and subcluster)
    # "" = source of data should be automatically determined (optional)
    if (len(missing_req_fields) == 0 and len(extra_fields) == 0):

        # Convert None values to "".
        set_empty(data_dict)
        # Confirm that certain keywords are lowercase.
        set_keywords(data_dict, keywords)
        # Determine which fields are not present, and add them to the dict.
        set_missing_keys(data_dict, optional_keys)
        # Certain fields should be lowercase.
        data_dict["type"] = data_dict["type"].lower()
        data_dict["description_field"] = data_dict["description_field"].lower()
        data_dict["eval_mode"] = data_dict["eval_mode"].lower()
        # For fields that are empty (""), set them with default values.
        set_dict_value(data_dict, "host_genus", "retrieve", "retain")
        set_dict_value(data_dict, "cluster", "retrieve", "retain")
        set_dict_value(data_dict, "subcluster", "retrieve", "retain")
        set_dict_value(data_dict, "annotation_author", "1", "retain")
        set_dict_value(data_dict, "retrieve_record", "1", "retain")
        set_dict_value(data_dict, "annotation_status", "draft", "final")

        # If status = 'draft', accession can remain "".
        if data_dict["annotation_status"] != "draft":
            set_dict_value(data_dict, "accession", "retrieve", "retain")
        return True
    else:
        print(f"The ticket {data_dict} is not formatted correctly.")
        return False



def parse_import_ticket_data(data_dict):
    """Converts import ticket data to a ImportTicket object.

    :param data_dict:
        A dictionary of data with the following keys:

            0. Import action type
            1. Primary PhageID
            2. Host
            3. Cluster
            4. Subcluster
            5. Status
            6. Annotation Author (int)
            7. Feature field
            8. Accession
            9. Retrieve Record (int)
            10. Eval mode

    :type data_dict: dict
    :returns: A pdm_utils ImportTicket object.
    :rtype: ImportTicket
    """
    ticket_attributes = constants.IMPORT_TABLE_STRUCTURE["valid_ticket"]
    other_attributes = data_dict.keys() - ticket_attributes

    tkt = ticket.ImportTicket()
    tkt.data_dict = data_dict
    for attr in ticket_attributes:
        attr_value = data_dict[attr]
        setattr(tkt, attr, attr_value)

    data_retrieve = set()
    data_retain = set()
    data_add = set()
    data_parse = set()

    other_attributes = list(other_attributes)
    x = 0
    while x < len(other_attributes):
        attr = other_attributes[x]
        attr_value = data_dict[attr]
        if attr_value == "retrieve":
            data_retrieve.add(attr)
        elif attr_value == "retain":
            data_retain.add(attr)
        elif attr_value == "parse":
            data_parse.add(attr)
        else:
            data_add.add(attr)
        x += 1
    tkt.data_retrieve = data_retrieve
    tkt.data_retain = data_retain
    tkt.data_parse = data_parse
    tkt.data_add = data_add

    return tkt


def set_empty(data_dict):
    """Convert None values to an empty string.

    :param data_dict: Dictionary of import ticket data.
    :type data_dict: dict
    """
    for key in data_dict.keys():
        if data_dict[key] is None:
            data_dict[key] = ""

def set_keywords(data_dict, keywords):
    """Convert specific values in a dictionary to lowercase.

    :param data_dict: Dictionary of import ticket data.
    :type data_dict: dict
    :param keywords:
        Set of valid keyword values that are handled differently
        than other values.
    :type keywords: set
    """
    for key in data_dict.keys():
        value = data_dict[key]
        if isinstance(value, str):
            if value.lower() in keywords:
                data_dict[key] = value.lower()


def set_missing_keys(data_dict, expected_keys):
    """Add a list of keys-values to a dictionary if it doesn't have those keys.

    :param data_dict: Dictionary of import ticket data.
    :type data_dict: dict
    :param expected_keys: Set of keys expected to be in the dictionary.
    :type expected_keys: set
    """
    missing_keys = expected_keys - data_dict.keys()
    for key in missing_keys:
        data_dict[key] = ""

def set_dict_value(data_dict, key, first, second):
    """Set the value for a specific key based on 'type' key-value.

    :param data_dict: Dictionary of import ticket data.
    :type data_dict: dict
    :param key: Dictionary key to change value of.
    :type key: str
    :param first: Value to assign to 'key' if 'type' == 'add'.
    :type first: str
    :param second: Value to assign to 'key' if 'type' != 'add'.
    :type second: str
    """
    if data_dict[key] == "":
        if data_dict["type"] == "add":
            data_dict[key] = first
        else:
            data_dict[key] = second

# TODO it probably makes more sense to pass the 'eval_mode' and
# 'description_field' data into modify_import_data() so that they get added
# to the data_dict. As a result, those values will be exported to the
# revised import table file when the script completes.
def construct_tickets(list_of_data_dict, eval_data_dict, description_field,
                      required_keys, optional_keys, keywords):
    """Construct pdm_utils ImportTickets from parsed data dictionaries.

    :param list_of_data_dict: List of import ticket data dictionaries.
    :type list_of_data_dict: list
    :param eval_data_dict: Dictionary of boolean evaluation flags.
    :type eval_data_dict: dict
    :param description_field:
        Default value to set ticket.description_field attribute if not
        present in the data dictionary.
    :type description_field: str
    :param required_keys: Set of keys required to be in the data dictionary.
    :type required_keys: set
    :param optional_keys:
        Set of optional keys that are not required to be in the data dictionary.
    :type optional_keys: set
    :param keywords:
        Set of valid keyword values that are handled differently
        than other values.
    :type keywords: set
    :returns: List of pdm_utils ImportTicket objects.
    :rtype: list
    """

    tkt_id = 0
    list_of_tickets = []
    for dict in list_of_data_dict:
        tkt_id += 1
        # Each ticket should contain a distinct eval_flag dictionary object.
        input_eval_mode = eval_data_dict["eval_mode"]
        input_eval_flag_dict = eval_data_dict["eval_flag_dict"].copy()

        result = modify_import_data(dict, required_keys, optional_keys, keywords)
        if result:
            tkt = parse_import_ticket_data(dict)
            tkt.id = tkt_id

            # Only set description_field and eval_mode from parameters
            # if they weren't set within the ticket.
            if tkt.description_field == "":
                tkt.description_field = description_field
                # Add description_field to data_dict so that it
                # can be saved to file.
                tkt.data_dict["description_field"] = description_field

            if tkt.eval_mode == "":
                tkt.eval_mode = input_eval_mode
                tkt.eval_flags = input_eval_flag_dict
                # Add eval_mode to data_dict so that it can be saved to file.
                tkt.data_dict["eval_mode"] = input_eval_mode

            else:
                if tkt.eval_mode in eval_modes.EVAL_MODES.keys():
                    tkt.eval_flags = eval_modes.get_eval_flag_dict(tkt.eval_mode)
            list_of_tickets.append(tkt)
        else:
            print("Unable to create ticket.")
    return list_of_tickets



def identify_duplicates(list_of_tickets, null_set=set()):
    """Compare all tickets to each other to identify ticket conflicts.

    Identifies if the same id, PhageID, and Accession
    is present in multiple tickets.

    :param list_of_tickets:
        A list of pdm_utils ImportTicket objects.
    :type list_of_tickets: list
    :param null_set:
        A set of values that may be expected to be duplicated, that
        should not throw errors.
    :type null_set: set
    :returns:
        tuple (tkt_id_dupes, phage_id_dupes)
        WHERE
        tkt_id_dupes(set) is a set of duplicate ticket ids.
        phage_id_dupes(set) is a set of duplicate PhageIDs.
    :rtype: tuple
    """
    tkt_id_list = []
    phage_id_list = []

    # Create separate lists to check each field for duplications.
    # Skip "none" values since they are expected to be duplicated.
    for tkt in list_of_tickets:
        if tkt.id not in null_set:
            tkt_id_list.append(tkt.id)
        if tkt.phage_id not in null_set:
            phage_id_list.append(tkt.phage_id)
    tkt_id_dupe_set = basic.identify_one_list_duplicates(tkt_id_list)
    phage_id_dupe_set = basic.identify_one_list_duplicates(phage_id_list)
    return (tkt_id_dupe_set, phage_id_dupe_set)


def get_genome(tkt, gnm_type=""):
    """Construct a pdm_utils Genome object from a pdm_utils ImportTicket object.

    :param tkt: A pdm_utils ImportTicket object.
    :type tkt: ImportTicket
    :param gnm_type: Identifier for the type of genome.
    :type gnm_type: str
    :returns: A pdm_utils Genome object.
    :rtype: Genome
    """
    gnm = genome.Genome()
    gnm.type = gnm_type
    gnm.set_id(value=tkt.phage_id)
    gnm.name = tkt.phage_id
    if "host_genus" in tkt.data_add:
        gnm.set_host_genus(tkt.data_dict["host_genus"])
    if "accession" in tkt.data_add:
        gnm.set_accession(tkt.data_dict["accession"])
    if "annotation_status" in tkt.data_add:
        gnm.annotation_status = tkt.data_dict["annotation_status"]
    if "cluster" in tkt.data_add:
        gnm.set_cluster(tkt.data_dict["cluster"])
    if "subcluster" in tkt.data_add:
        gnm.set_subcluster(tkt.data_dict["subcluster"])
    if "annotation_author" in tkt.data_add:
        gnm.set_annotation_author(tkt.data_dict["annotation_author"])
    if "retrieve_record" in tkt.data_add:
        gnm.set_retrieve_record(tkt.data_dict["retrieve_record"])
    return gnm
