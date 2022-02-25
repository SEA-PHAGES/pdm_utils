"""Configuration file definition and parsing."""

import configparser
import sys

from pdm_utils.functions import basic


def default_sections_keys():
    sections = {"mysql": {"user", "password"},
                "mysqldump": {"user", "password"},
                "ncbi": {"api_key", "email", "tool"},
                "upload_server": {"host", "dest", "user", "password"},
                "download_server": {"url"},
                "emailer": {"username", "password"}}
    return sections


def setup_section(keys, value):
    section = {}
    for key in keys:
        section[key] = value
    return section


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
    if file is not None:
        parser = parse_config(file, parser)
    return parser


def write_config(parser, filepath):
    """Write a ConfigParser to file."""
    with filepath.open("w") as fh:
        parser.write(fh)


def create_empty_config_file(dir, file, null_value):
    """Create an empty config file with all available settings."""
    output_path = basic.set_path(dir, kind="dir", expect=True)
    config_path = basic.make_new_file(output_path, file, "txt", attempt=50)
    if config_path is None:
        print("Unable to create config file. File already exists.")
    else:
        parser = default_parser(null_value)
        write_config(parser, config_path)
