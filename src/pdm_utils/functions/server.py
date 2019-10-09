"""Misc. functions to utilize webfaction server."""

import getpass
from pdm_utils.constants import constants
import paramiko

# Add log file for paramiko because it throws an error otherwise.
# Soft requirement for compliance with Paramiko standards
paramiko.util.log_to_file("/tmp/paramiko.log")


# TODO unittest (but manually tested).
def get_credentials():
    """Get server credentials."""
    username = getpass.getpass(prompt="Server username: ")
    password = getpass.getpass(prompt="Server password: ")
    return (username, password)

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
def setup_connection(transport, attempts=3):
    """Get credentials and setup connection to the server."""
    sftp = None
    while attempts > 0 and sftp is None:
        user, pwd = get_credentials()
        try:
            transport.connect(username=user,password=pwd)
            sftp = paramiko.SFTPClient.from_transport(transport)
        except:
            print("Unable to connect to server."
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
        print("Problems encountered with the"
              "specified local or destination directories.")
        result = False
    return result
