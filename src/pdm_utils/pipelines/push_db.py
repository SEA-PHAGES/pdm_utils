"""Pipeline to push data to the server."""
import argparse
import pathlib
import sys

from pdm_utils.functions import basic
from pdm_utils.constants import constants
from pdm_utils.functions import configfile
from pdm_utils.functions import server


# TODO unittest.
def main(unparsed_args_list):
    """Run the push_db pipeline."""
    args = parse_args(unparsed_args_list)
    local_dir = args.directory
    local_file = args.file

    # Create config object with data obtained from file and/or defaults.
    # For server host and dir, give priority to command line over config file.
    # So that overriding "default" option is easier
    config = configfile.build_complete_config(args.config_file)
    server_host = config["upload_server"]["host"]
    remote_dir = config["upload_server"]["dest"]
    user = config["upload_server"]["user"]
    pwd = config["upload_server"]["password"]

    if server_host is None or args.server_host is not None:
        server_host = args.server_host
    if remote_dir is None or args.remote_directory is not None:
        remote_dir = args.remote_directory

    if server_host is None or remote_dir is None:
        print("No host and/or remote directory provided. "
              "Unable to upload file(s).")
        sys.exit(1)

    file_list = get_files(local_dir, local_file, set([".DS_Store"]))

    status = True
    if len(file_list) == 0:
        print("There are no files to upload.")
        status = False

    if status is True:
        server.set_log_file(str(args.log_file))
        transport = server.get_transport(server_host)
        if transport is None:
            status = False

    if status is True:
        sftp = server.setup_sftp_conn(transport, user, pwd, attempts=1)
        if sftp is None:
            status = False

    if status is True:
        success, fail = upload(sftp, remote_dir, file_list)
        sftp.close()
        transport.close()

        if len(fail) > 0:
            print("The following files were not uploaded:")
            for file in fail:
                print(file)


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for uploading to the server."""

    push_db_help = "Pipeline to upload files to a server."
    directory_help = "Path to the folder containing files for upload."
    file_help = "Path to the file for upload."
    log_file_help = "Path to the file to log paramiko output."
    server_host_help = "Server host name."
    remote_directory_help = "Server directory to upload files."
    config_file_help = "Path to the file containing user-specific login data."

    parser = argparse.ArgumentParser(description=push_db_help)
    parser.add_argument("-d", "--directory", type=pathlib.Path,
                        help=directory_help)
    parser.add_argument("-f", "--file", type=pathlib.Path,
                        help=file_help)
    parser.add_argument("-l", "--log_file", type=pathlib.Path,
                        default=pathlib.Path("/tmp/paramiko.log"),
                        help=log_file_help)
    parser.add_argument("-s", "--server_host", type=str,
                        default=None, help=server_host_help)
    parser.add_argument("-rd", "--remote_directory", type=str,
                        default=None, help=remote_directory_help)
    parser.add_argument("-c", "--config_file", type=pathlib.Path,
                        help=config_file_help, default=None)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args


# TODO test.
def get_files(directory, file, ignore_set):
    """Get the list of file(s) that need to be uploaded."""
    file_list = []

    if directory is not None:
        directory = basic.set_path(directory, kind="dir", expect=True)
        folder_files = basic.identify_contents(directory, kind="file",
                                               ignore_set=ignore_set)
        file_list.extend(folder_files)

    if file is not None:
        file = basic.set_path(file, kind="file", expect=True)
        file_list.append(file)

    return file_list


# TODO test.
def upload(sftp, remote_dir, file_list):
    """Upload file(s)."""
    success = []
    fail = []
    for local_filepath in file_list:
        print(f"Uploading {local_filepath.name}...")
        remote_filepath = pathlib.Path(remote_dir, local_filepath.name)
        result = server.upload_file(sftp, str(local_filepath),
                                    str(remote_filepath))
        if result:
            success.append(local_filepath.name)
        else:
            fail.append(local_filepath.name)
    return success, fail
