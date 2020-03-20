"""Eval mode functions and dictionaries."""

from pdm_utils.functions import basic

# Eval_flag definitions.
EVAL_FLAGS = {
    # Options that impact how data is processed but are not utilized
    # for evaluation of specific parts of a flat file:
    "check_replace": "Should unexpected genome replacements be reported? ",
    "import_locus_tag": "Should CDS feature locus_tags be imported? ",

    # Options that are utilized during the evaluation stage:
    "check_locus_tag": "Should the structure of CDS feature locus_tags be evaluated? ",
    "check_description_tally": "Should the number of descriptions be evaluated? ",
    "check_description_field": "Should CDS descriptions in unexpected fields be reported? ",
    "check_description": "Should unexpected CDS descriptions be reported? ",
    "check_trna": "Should tRNA features be evaluated? ",
    "check_id_typo": "Should genome ID typos be reported? ",
    "check_host_typo": "Should host typos be reported? ",
    "check_author": "Should unexpected authors be reported? ",
    "check_gene": "Should the CDS 'gene' qualifier be evaluated? ",
    "check_seq": "Should the nucleotide sequence be evaluated? ",
    "check_coords": "Should feature duplication be evaluated? "
    }

# Eval mode definitions.
EVAL_MODES = {
    "draft": ("Relaxed evaluations for automatically generated draft "
              "genome annotations since data have not been manually reviewed."),
    "final": ("Stringent evaluations for manual genome annotations."),
    "auto": ("Relaxed evaluations for genome annotations "
             "automatically retrieved from GenBank for routine replacement."),
    "misc": ("Relaxed evaluations for genome annotations "
             "that have been generated from external sources."),
    "custom": ("User-defined evaluations for customized import.")
    }

def get_eval_flag_dict(eval_mode):
    """Get a dictionary of evaluation flags.

    :param eval_mode:
        Valid evaluation mode (base, draft, final, auto, misc, custom)
    :type eval_mode: str
    :returns: Dictionary of boolean values.
    :rtype: dict
    """

    # Base dictionary with all flags set to True.
    dict = {}
    for key in EVAL_FLAGS:
        dict[key] = True

    # Auto-annotations.
    if eval_mode == "draft":
        dict["check_locus_tag"] = False
        dict["check_trna"] = False
        dict["import_locus_tag"] = False
        dict["check_id_typo"] = False
        dict["check_host_typo"] = False
        dict["check_author"] = False
        dict["check_description"] = False
        dict["check_coords"] = False

    # Manual annotations.
    elif eval_mode == "final":
        dict["import_locus_tag"] = False

    # SEA-PHAGES GenBank records.
    elif eval_mode == "auto":
        dict["check_locus_tag"] = False
        dict["check_description_field"] = False
        dict["check_replace"] = False
        dict["check_trna"] = False
        dict["check_id_typo"] = False
        dict["check_host_typo"] = False
        dict["check_author"] = False
        dict["check_description"] = False
        dict["check_description_tally"] = False
        dict["check_gene"] = False
        dict["check_coords"] = False

    # Non-SEA-PHAGES GenBank records.
    elif eval_mode == "misc":
        dict["check_locus_tag"] = False
        # TODO below should probably be True, but it causes problems
        # when checking the current genome, GNM2_001, since these are not 'draft'
        # genomes.
        dict["check_replace"] = False
        dict["check_trna"] = False
        dict["check_id_typo"] = False
        dict["check_host_typo"] = False
        dict["check_author"] = False
        dict["check_description"] = False
        dict["check_description_tally"] = False
        dict["check_gene"] = False

    # Custom QC settings. User can select the settings, so it is initialized as
    # a copy of the base eval_mode. The user can provide the
    # customized combination of options.
    elif eval_mode == "custom":
        for key in dict.keys():
            prompt = f"Eval_flag: {key}. {EVAL_FLAGS[key]}"
            response = basic.ask_yes_no(prompt=prompt, response_attempt=3)
            if response is None:
                print("The default setting for this eval_flag will be used.")
            else:
                dict[key] = response

    elif eval_mode == "base":
        pass
    else:
        print("A valid eval_mode has not been selected.")
    return dict
