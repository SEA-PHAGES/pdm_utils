"""Pipeline to push data to the server."""
import argparse
import pathlib
from pdm_utils.constants import constants
from pdm_utils.functions import basic
from pdm_utils.functions import server


# TODO unittest.
def main(unparsed_args_list):
    """Run the push_db pipeline."""
    args = parse_args(unparsed_args_list)

    file_list = []
    if args.directory is not None:
        args.directory = basic.set_path(args.directory,
                                              kind="dir", expect=True)
        folder_files = basic.identify_contents(args.directory, kind="file",
                                               ignore_set=set([".DS_Store"]))
        file_list.extend(folder_files)
    if args.file is not None:
        args.file = basic.set_path(args.file, kind="file", expect=True)
        file_list.append(args.file)

    status = True
    if len(file_list) == 0:
        print("There are no files to upload.")
        status = False

    if status == True:
        server.set_log_file(str(args.log_file))
        transport = server.get_transport(constants.DB_HOST)
        if transport is None:
            status = False

    if status == True:
        sftp = server.setup_sftp_conn(transport, attempts=3)
        if sftp is None:
            status = False

    success = []
    fail = []
    if status == True:
        for local_filepath in file_list:
            print(f"Uploading {local_filepath.name}...")
            remote_filepath = pathlib.Path(constants.DB_HOST_DIR,
                                           local_filepath.name)
            result = server.upload_file(sftp, str(local_filepath),
                                        str(remote_filepath))
            if result:
                success.append(local_filepath.name)
            else:
                fail.append(local_filepath.name)
        sftp.close()
        transport.close()

    if len(fail) > 0:
        print("The following files were not uploaded:")
        for file in fail:
            print(file)


# TODO unittest.
def parse_args(unparsed_args_list):
    """Verify the correct arguments are selected for uploading to the server."""
    PUSH_DB_HELP = ("Pipeline to upload a new version of "
                   "a MySQL database to the server.")
    DIRECTORY_HELP = ("Path to the folder containing files for upload.")
    FILE_HELP = ("Path to the file for upload.")
    LOG_FILE_HELP = ("Path to the file to log paramiko output.")
    parser = argparse.ArgumentParser(description=PUSH_DB_HELP)
    parser.add_argument("-d", "--directory", type=pathlib.Path,
        help=DIRECTORY_HELP)
    parser.add_argument("-f", "--file", type=pathlib.Path,
        help=FILE_HELP)
    parser.add_argument("-l", "--log_file", type=pathlib.Path,
        default=pathlib.Path("/tmp/paramiko.log"), help=LOG_FILE_HELP)

    # Assumed command line arg structure:
    # python3 -m pdm_utils.run <pipeline> <additional args...>
    # sys.argv:      [0]            [1]         [2...]
    args = parser.parse_args(unparsed_args_list[2:])
    return args
