"""Configuration file definition and parsing."""

import configparser
import sys

from pdm_utils.functions import basic

CONFIG_DICT = {"mysql": {"user", "password"},
               "mysqldump": {"user", "password"},
               "ncbi": {"ncbi_api_key", "ncbi_email", "ncbi_tool"}
               }

# TODO test.
def subconfig(keys):
    dict = {}
    # Default value for everything is "", since configparser requires strings.
    for key in keys:
        dict[key] = ""
    return dict

# TODO test.
def default_config():
    """Constructs complete config with empty values."""
    null_config = configparser.ConfigParser()
    for key in CONFIG_DICT.keys():
        null_config[key] = subconfig(CONFIG_DICT[key])
    return null_config

# TODO test.
def parse_config(file):
    """Get parameters from config file."""
    filepath = basic.set_path(file, kind="file", expect=True)
    config = configparser.ConfigParser()
    try:
        config.read(filepath)
    except:
        print("Unable to parse config file")
        sys.exit(1)
    else:
        return config

# TODO test.
def build_complete_config(file):
    "Buid a complete config object by merging user-supplied and default config."
    # The config file may not have all possible values.
    # So parse and get all data in the config file and merge it to a default
    # null config object.
    input_config = parse_config(file)
    full_config = default_config()
    for section in full_config.keys():
        if section in input_config.keys():
            for key in full_config[section].keys():
                try:
                    full_config[section][key] = input_config[section][key]
                except:
                    pass
    return full_config

# TODO test.
def reformat_data(config_dict, old, new):
    """Convert string value to another value."""
    dict = {}
    for key in config_dict.keys():
        if config_dict[key] == old:
            dict[key] = new
        else:
            dict[key] = config_dict[key]
    return dict
