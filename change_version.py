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

VERSION1 = '__version__ = "{}"\n'
VERSION2 = "version = '{}'\n"
VERSION3 = "release = '{}'\n"
VERSION4 = '    version="{}",\n'
VERSION_FILES = {
    "init_file":{
        "path":INIT_PATH,
        "lines":{
            13:{VERSION1[:13]:VERSION1}
        }
    },
    "conf_file":{
        "path":CONF_PATH,
        "lines":{
            9:{VERSION2[:9]:VERSION2,
               VERSION3[:9]:VERSION3}
        }
    },
    "setup_file":{
        "path":SETUP_PATH,
        "lines":{
            12: {VERSION4[:12]:VERSION4}
        }
    }
}


def main(unparsed_args):
    args = parse_args(unparsed_args)
    versions = get_all_versions(VERSION_FILES)
    old_version = compare_versions(versions)
    old_components = split_components(old_version)
    if args.component == "check":
        print("\n\nNo changes made to version.")
    else:
        new_components = change_components(old_components, args.component,
                                           args.increment)
        new_version = join_components(new_components)
        print("\n\nChanging version in all locations...")
        print(f"\nNew version: {new_version}")
        print("\nChanged lines:\n")
        for file_name in VERSION_FILES.keys():
            file_data = VERSION_FILES[file_name]
            write_out(file_data, new_version)
        print("\n\nVersion updated.")


def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for changing the version."""
    CHANGE_VERSION_HELP = \
        ("Script to change the version number in all "
        "locations in the repository at once. "
        "It confirms that each version line can be identified, "
        "each version line is structured as expected, "
        "and that all versions in all locations are identical. "
        "If selected, it then will change the version in each location.")
    COMPONENT_HELP = \
        ("The version component to change. "
         "'Check' will check all version data without making any changes.")
    INCREMENT_HELP = ("Amount that version component should be incremented.")
    parser = argparse.ArgumentParser(description=CHANGE_VERSION_HELP)
    parser.add_argument("component", type=str,
        choices=["major", "minor", "micro", "check"], help=COMPONENT_HELP)
    parser.add_argument("-i", "--increment", type=int,
        default=1, help=INCREMENT_HELP)

    # Assumed command line arg structure:
    # python3 ./control_version.py component -v 1
    # sys.argv:      [0]            [1]     [2] [3]
    args = parser.parse_args(unparsed_args_list[1:])
    return args


def get_all_versions(file_dict):
    all_versions = []
    for file_name in file_dict.keys():
        file_data = file_dict[file_name]
        file_versions = get_file_versions(file_data)
        all_versions.extend(file_versions)
    return all_versions


def get_file_versions(file_data):
    file_versions = []
    lines_found = set()
    file_path = file_data["path"]
    print(f"\n\nChecking {file_path.name} version data...")
    with file_path.open("r") as file_handle:
        data_list = file_handle.readlines()

    expected = get_expected_lines(file_data["lines"])
    actual = 0
    for string in data_list:
        for offset in file_data["lines"]:
            lines = file_data["lines"][offset]
            if string[:offset] in lines.keys():
                lines_found.add(string[:offset])
                actual += 1
                version = string[offset:].strip()
                version = version.replace('"', '')
                version = version.replace("'", '')
                version = version.replace(",", '')
                file_versions.append(version)

    if expected - actual != 0:
        print("Error finding version lines.")
        print(f"\n{actual} line(s) found:")
        for line in lines_found:
            print(line)

        print(f"\n{expected} line(s) expected:")
        for offset in file_data["lines"]:
            lines = file_data["lines"][offset]
            for line in lines.keys():
                print(line)
        sys.exit(1)
    else:
        return file_versions


def get_expected_lines(line_data):
    lines = 0
    for offset in line_data:
        lines += len(line_data[offset].keys())
    return lines


def compare_versions(versions):
    version_set = set(versions)
    if len(version_set) != 1:
        print("\nMultiple versions identified:")
        for version in version_set:
            print("'" + version + "'")
        sys.exit(1)
    else:
        version = list(version_set)[0]
        print("\nAll versions are identical:")
        print("'" + version + "'")
        return version


def split_components(version):
    version_data = version.split(".")
    major = int(version_data[0])
    minor = int(version_data[1])
    micro = int(version_data[2])
    return (major, minor, micro)


def change_components(components, change_type, increment):
    major = components[0]
    minor = components[1]
    micro = components[2]
    if change_type == "major":
        new_major = update_value(major, increment)
        new_minor = update_value(minor, -1 * minor)
        new_micro = update_value(micro, -1 * micro)
    elif change_type == "minor":
        new_major = major
        new_minor = update_value(minor, increment)
        new_micro = update_value(micro, -1 * micro)
    else:
        new_major = major
        new_minor = minor
        new_micro = update_value(micro, increment)

    return (new_major, new_minor, new_micro)


def update_value(current, delta):
    new = current + delta
    return new


def join_components(components):
    major = components[0]
    minor = components[1]
    micro = components[2]
    version = ".".join([str(major), str(minor), str(micro)])
    return version


def write_out(file_data, new_version):
    file_path = file_data["path"]
    with file_path.open("r") as file_handle:
        data_list = file_handle.readlines()
    index = 0
    while index < len(data_list):
        string = data_list[index]
        for offset in file_data["lines"]:
            lines = file_data["lines"][offset]
            if string[:offset] in lines.keys():
                for key in lines.keys():
                    if string[:offset] == key:
                        string = lines[key].format(new_version)
                        print(string)
        data_list[index] = string
        index += 1
    with file_path.open("w") as file_handle:
         file_handle.writelines(data_list)


if __name__ == "__main__":
    main(sys.argv)
