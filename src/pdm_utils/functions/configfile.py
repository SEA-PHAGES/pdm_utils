"""Configuration file definition and parsing."""

import configparser
import sys

from pdm_utils.functions import basic

def default_sections_keys():
    dict = {"mysql": {"user", "password"},
            "mysqldump": {"user", "password"},
            "ncbi": {"api_key", "email", "tool"},
            "upload_server": {"host", "dest", "user", "password"}
            }
    return dict

def setup_section(keys, value):
    dict = {}
    for key in keys:
        dict[key] = value
    return dict

def default_parser(null_value):
    """Constructs complete config with empty values."""
    # Need to allow no value if the null value is None.
    null_parser = configparser.ConfigParser(allow_no_value=True)
    all_configs = default_sections_keys()
    for key in all_configs.keys():
        null_parser[key] = setup_section(all_configs[key], null_value)
    return null_parser

def parse_config(file, parser=None):
    """Get parameters from config file."""
    filepath = basic.set_path(file, kind="file", expect=True)
    if parser is None:
        parser = configparser.ConfigParser()

    try:
        parser.read(filepath)
    except:
        print("Unable to parse config file")
        sys.exit(1)
    else:
        return parser

def build_complete_config(file):
    "Buid a complete config object by merging user-supplied and default config."
    parser = default_parser(None)
    parser = parse_config(file, parser)
    return parser
