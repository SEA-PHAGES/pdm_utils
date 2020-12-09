"""
Pipeline to push files to a server using SFTP.
"""

import argparse
import getpass
import pathlib
import sys

import paramiko

from pdm_utils.functions import basic
from pdm_utils.functions import configfile


# TODO unittest.
def parse_args(unparsed_args):
    """
    Verify the correct arguments are selected for uploading to the server.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--directory",
                        type=pathlib.Path,
                        help="path to a folder whose contents will be "
                             "uploaded to the server")
    parser.add_argument("-f", "--file",
                        type=pathlib.Path,
                        help="path to a file to be uploaded to the server")
    parser.add_argument("-l", "--log_file",
                        type=pathlib.Path,
                        default=pathlib.Path("/tmp/paramiko.log"),
                        help="path to a file where paramiko output should "
                             "be logged")
    parser.add_argument("-s", "--server_host",
                        type=str,
                        default=None,
                        help="server hostname")
    parser.add_argument("-k", "--key-file",
                        type=pathlib.Path,
                        help="path to a key-file to use for server login")
    parser.add_argument("-rd", "--remote_directory",
                        type=pathlib.Path,
                        default=None,
                        help="path to server directory where uploaded file(s) "
                             "should go")
    parser.add_argument("-c", "--config_file",
                        type=pathlib.Path,
                        default=None,
                        help="path to local file containing server login data")

    args = parser.parse_args(unparsed_args)
    return args


# TODO test.
def get_files(directory, file, ignore):
    """
    Get the list of file(s) that need to be uploaded.

    :param directory: (optional) directory containing files for upload
    :type: directory: pathlib.Path
    :param file: (optional) file to upload
    :type file: pathlib.Path
    :param ignore: file(s) to ignore during upload process
    :type ignore: set
    :return: file_list
    """
    file_list = []

    if directory is not None:
        directory = basic.set_path(directory, kind="dir", expect=True)
        folder_files = basic.identify_contents(directory, kind="file",
                                               ignore_set=ignore)
        file_list.extend(folder_files)

    if file is not None:
        file = basic.set_path(file, kind="file", expect=True)
        file_list.append(file)

    return file_list


def upload(sftp_client, destination, files):
    """
    Try to upload the file(s).

    :param sftp_client: an open SFTPClient
    :type sftp_client: paramiko.SFTPClient
    :param destination: remote file directory to upload to
    :type destination: pathlib.Path
    :param files: the file(s) to upload
    :type files: list of pathlib.Path
    :return: successes, failures
    """
    successes, failures = list(), list()

    for local_file in files:
        name = local_file.name
        print(f"Uploading {name}...")
        remote_file = destination.joinpath(name)
        try:
            sftp_client.put(str(local_file), str(remote_file))
            successes.append(name)
        except OSError:
            failures.append(name)

    return successes, failures


# TODO unittest.
def main(unparsed_args):
    """
    Driver function for the push pipeline.

    :param unparsed_args: the command-line arguments given to this
                          pipeline's caller (likely pdm_utils.__main__)
    :type unparsed_args: list
    """
    # Parse the command line args with argparse
    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parse_args(unparsed_args[2:])

    # Parse config file if one was given
    config = configfile.build_complete_config(args.config_file)
    server_host = config["upload_server"]["host"]
    remote_dir = config["upload_server"]["dest"]
    user = config["upload_server"]["user"]
    password = config["upload_server"]["password"]

    # Command line hostname and destination override config file defaults
    if args.server_host is not None:
        server_host = args.server_host
    if args.remote_directory is not None:
        remote_dir = args.remote_directory
    key_file = args.key_file

    # Can't upload files to unknown host or destination
    if server_host is None or remote_dir is None:
        print("No hostname and/or remote directory provided. Unable "
              "to upload file(s).")
        sys.exit(1)

    # Get the list of files to upload
    file_list = get_files(args.directory, args.file, ignore={".DS_Store"})

    if len(file_list) > 0:
        if args.log_file:
            paramiko.util.log_to_file(args.log_file)

        # Keep track of failed uploads
        failures = list()

        # Setup paramiko connection to server
        with paramiko.Transport(server_host) as transport:
            # Context manager so Transport will be closed automatically
            if user and key_file:
                # Username and key-file is preferred authentication method
                try:
                    key = paramiko.RSAKey.from_private_key_file(key_file)
                except paramiko.ssh_exception.SSHException:
                    print(f"'{key_file}' is not a valid RSA key file")
                    sys.exit(1)
                try:
                    transport.connect(username=user, pkey=key)
                except paramiko.ssh_exception.AuthenticationException:
                    print(f"Authentication failed with user '{user}' and "
                          f"key-file '{key_file}'")
                    sys.exit(1)
            elif user and password:
                # Username and password from config is next preferred method
                try:
                    transport.connect(username=user, password=password)
                except paramiko.ssh_exception.AuthenticationException:
                    p = f"{password[0]}{'*' * (len(password)-2)}{password[-1]}"
                    print(f"Authentication failed with user '{user}' and "
                          f"password '{p}'")
                    sys.exit(1)
            else:
                # Finally, get username and password from command line
                user = getpass.getpass(f"Enter username for {server_host}: ")
                password = getpass.getpass(f"Enter password for "
                                           f"{user}@{server_host}: ")
                try:
                    transport.connect(username=user, password=password)
                except paramiko.ssh_exception.AuthenticationException:
                    p = f"{password[0]}{'*' * (len(password)-2)}{password[-1]}"
                    print(f"Authentication failed with user '{user}' and "
                          f"password '{p}'")
                    sys.exit(1)
            with paramiko.SFTPClient.from_transport(transport) as sftp_client:
                # Context manager so SFTPClient will be closed automatically
                for local_file in file_list:
                    # sftp_client.put requires remote filename not remote dir
                    remote_file = remote_dir.joinpath(local_file.name)
                    try:
                        print(f"Uploading {str(local_file)}...")
                        sftp_client.put(str(local_file), str(remote_file))
                    except OSError:
                        failures.append(local_file)
        for file in failures:
            print(f"Could not upload {str(file)}")
    else:
        print("No files to upload")
