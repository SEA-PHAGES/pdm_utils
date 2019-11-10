"""Run mode functions and dictionaries."""

from pdm_utils.functions import basic

# Eval_flag definitions.
EVAL_FLAGS = {
    # Options that impact how data is processed but are not utilized
    # for evaluation of specific parts of a flat file:
    "check_replace": "Should unexpected genome replacements be reported? ",
    "import_locus_tag": "Should CDS feature locus_tags be imported? ",

    # Options that are utilized during the evaluation stage:
    "check_locus_tag": "Should the structure of CDS feature locus_tags be checked? ",
    "check_description_field": "Should CDS descriptions in unexpected fields be reported? ",
    "check_description": "Should unexpected CDS descriptions be reported? ",
    "check_trna": "Should tRNA features be evaluated? ",
    "check_id_typo": "Should genome ID typos be reported? ",
    "check_host_typo": "Should host typos be reported? ",
    "check_author": "Should unexpected authors be reported? ",
    "check_gene": "Should the CDS 'gene' qualifier be evaluated? ",
    "check_seq": "Should the nucleotide sequence be evaluated? "
    }

# Run mode definitions.
RUN_MODES = {
    "pecaan": ("Relaxed evaluations for draft genome annotations "
               "retrieved from PECAAN since some data has not yet been "
               "manually reviewed (such as locus_tags)."),
    "phagesdb": ("Most stringent evaluations for final genome annotations "
                 "retrieved from PhagesDB since this data will be "
                 "submitted to GenBank."),
    "sea_auto": ("Relaxed evaluations for genome annotations generated "
                 "through SEA-PHAGES but retrieved from GenBank since "
                 "some of this data can no longer be modified."),
    "misc": ("Very relaxed evaluations for genome annotations "
             "not generated through SEA-PHAGES, since the data cannot be "
             "modified."),
    "custom": ("User-defined evaluations for customized import.")
    }

def get_eval_flag_dict(run_mode):
    """."""
    # Base dictionary with all flags set to True.
    dict = {}
    for key in EVAL_FLAGS:
        dict[key] = True

    # Auto-annotations.
    if run_mode == "pecaan":
        dict["check_locus_tag"] = False
        dict["check_trna"] = False
        dict["import_locus_tag"] = False
        dict["check_id_typo"] = False
        dict["check_host_typo"] = False
        dict["check_author"] = False
        dict["check_description"] = False

    # Manual annotations.
    elif run_mode == "phagesdb":
        dict["import_locus_tag"] = False

    # SEA-PHAGES GenBank records.
    elif run_mode == "sea_auto":
        dict["check_locus_tag"] = False
        dict["check_description_field"] = False
        dict["check_replace"] = False
        dict["check_trna"] = False
        dict["check_id_typo"] = False
        dict["check_host_typo"] = False
        dict["check_author"] = False
        dict["check_description"] = False
        dict["check_gene"] = False

    # Non-SEA-PHAGES GenBank records.
    elif run_mode == "misc":
        dict["check_locus_tag"] = False
        dict["check_replace"] = False
        dict["check_trna"] = False
        dict["check_id_typo"] = False
        dict["check_host_typo"] = False
        dict["check_author"] = False
        dict["check_description"] = False
        dict["check_gene"] = False

    # Custom QC settings. User can select the settings, so it is initialized as
    # a copy of the base run_mode. The user can provide the
    # customized combination of options.
    elif run_mode == "custom":
        for key in dict.keys():
            prompt = f"Eval_flag: {key}. {EVAL_FLAGS[key]}"
            response = basic.ask_yes_no(prompt=prompt, response_attempt=3)
            if response is None:
                print("The default setting for this eval_flag will be used.")
            else:
                dict[key] = response
    else:
        print("A valid run_mode has not been selected.")
    return dict
