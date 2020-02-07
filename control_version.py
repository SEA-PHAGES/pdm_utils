"""Script to change the package version number in all locations at once."""
import argparse
import pathlib
import sys

# Path to all files containing version data.
PACKAGE_PATH = pathlib.Path(__file__)
PACKAGE_DIR = PACKAGE_PATH.parent
INIT_PATH = pathlib.Path(PACKAGE_DIR, "src/pdm_utils/__init__.py")
SETUP_PATH = pathlib.Path(PACKAGE_DIR, "setup.py")
CONF_PATH = pathlib.Path(PACKAGE_DIR, "docs/source/conf.py")

def update_value(current, delta):
    new = current + delta
    return new

def main(unparsed_args):
    args = parse_args(unparsed_args)
    version = get_current_version(INIT_PATH)
    major, minor, micro = split_components(version)
    if args.component == "major":
        new_major = update_value(major, args.value)
        new_minor = update_value(minor, -1 * minor)
        new_micro = update_value(micro, -1 * micro)
    elif args.component == "minor":
        new_major = major
        new_minor = update_value(minor, args.value)
        new_micro = update_value(micro, -1 * micro)
    else:
        new_major = major
        new_minor = minor
        new_micro = update_value(micro, args.value)

    new_version = join_components(new_major, new_minor, new_micro)
    write_out(INIT_PATH, new_version)
    write_out(SETUP_PATH, new_version)
    write_out(CONF_PATH, new_version)
    print(f"Original version: {version}")
    print(f"New version: {new_version}")

def get_current_version(file_path):
    with file_path.open("r") as file_handle:
        data_list = file_handle.readlines()
        for string in data_list:
            if string.startswith("__version__ ="):
                string_split = string.split("=")
                version = string_split[1].strip()
                version = version.replace('"', '')
    return version

def split_components(version):
    version_data = version.split(".")
    major = int(version_data[0])
    minor = int(version_data[1])
    micro = int(version_data[2])
    return major, minor, micro

def join_components(major, minor, micro):
    version = ".".join([str(major), str(minor), str(micro)])
    return version

def write_out(file_path, new_version):
    with file_path.open("r") as file_handle:
        data_list = file_handle.readlines()
    index = 0
    while index < len(data_list):
        string = data_list[index]
        if file_path.stem == "__init__":
            if string.startswith("__version__ ="):
                string = f'__version__ = "{new_version}"\n'
        if file_path.stem == "conf":
            if string.startswith("version ="):
                string = f"version = '{new_version}'\n"
                # input("pause")
            elif string.startswith("release ="):
                string = f"release = '{new_version}'\n"
                # input("pause2")
            else:
                pass
        if file_path.stem == "setup":
            if string.startswith("    version="):
                string = f'    version="{new_version}",\n'
        data_list[index] = string
        index += 1
    with file_path.open("w") as file_handle:
         file_handle.writelines(data_list)

def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for changing the version."""
    CONTROL_VERSION_HELP = ("Script to change the version number in all locations.")
    COMPONENT_HELP = ("Part of the version to change.")
    VALUE_HELP = ("Amt to change the component.")
    parser = argparse.ArgumentParser(description=CONTROL_VERSION_HELP)
    parser.add_argument("component", type=str,
        choices=["major", "minor", "micro"], help=COMPONENT_HELP)
    parser.add_argument("-v", "--value", type=int,
        default=1, help=VALUE_HELP)

    # Assumed command line arg structure:
    # python3 ./control_version.py component -v 1
    # sys.argv:      [0]            [1]     [2] [3]
    args = parser.parse_args(unparsed_args_list[1:])
    return args

if __name__ == "__main__":
    main(sys.argv)
