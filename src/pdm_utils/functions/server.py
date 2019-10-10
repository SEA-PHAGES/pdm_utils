"""Misc. functions to utilizes server."""

import getpass
import paramiko
from pdm_utils.functions import basic


# TODO unittest (but manually tested).
def set_log_file(filepath):
    """Set the filepath used to stored the Paramiko output.

    This is a soft requirement for compliance with Paramiko standards.
    If it is not set, paramiko throws an error.
    """
    paramiko.util.log_to_file(filepath)


# TODO unittest (but manually tested).
def get_transport(host):
    """Create paramiko Transport with the server name."""
    # If the host name is not valid, paramiko is unable to find the host
    # and throws an error.
    try:
        transport = paramiko.Transport(host)
    except:
        transport = None
        print("Unable to find server: " + host)
    return transport


# TODO unittest (but manually tested).
def setup_sftp_conn(transport, attempts=1):
    """Get credentials and setup connection to the server."""
    sftp = None
    # Note: the way paramiko manages the connectionn, this loop
    # doesn't seem to work with multiple attempts.
    # There may be a transport attribute that needs to be reset or
    # a new transport object may need to be generated after each failed
    # attempt.
    while (attempts > 0 and sftp is None):
        user, pwd = basic.get_user_pwd(user_prompt="Server username: ",
                                       pwd_prompt="Server password: ")
        try:
            transport.connect(username=user, password=pwd)
            sftp = paramiko.SFTPClient.from_transport(transport)
        except:
            print("Unable to connect to server. "
                  "Incorrect username and password")
            attempts -= 1
    return sftp

# TODO unittest (but manually tested).
def upload_file(sftp, local_filepath, remote_filepath):
    """Upload a file to the server."""
    try:
        sftp.put(local_filepath, remote_filepath)
        result = True
    except:
        print("Problems encountered with the "
              "specified local or destination directories.")
        result = False
    return result
