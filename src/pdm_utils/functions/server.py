"""Misc. functions to utilizes server."""

import getpass
import paramiko
from pdm_utils.functions import basic


# TODO unittest (but manually tested).
def set_log_file(filepath):
    """Set the filepath used to stored the Paramiko output.

    This is a soft requirement for compliance with Paramiko standards.
    If it is not set, paramiko throws an error.

    :param filepath: Path to file to log Paramiko results.
    :type filepath: Path
    """
    paramiko.util.log_to_file(filepath)


# TODO unittest (but manually tested).
def get_transport(host):
    """Create paramiko Transport with the server name.

    :param host: Server to connect to.
    :type host: str
    :returns:
        Paramiko Transport object. If the server is not available,
        None is returned.
    :rtype: Transport
    """
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
    """Get credentials and setup connection to the server.

    :param transport: Paramiko Transport object directed towards a valid server.
    :type transport: Transport
    :param attempts: Number of attempts to connect to the server.
    :type attempts: int
    :returns:
        Paramiko SFTPClient connection. If no connection can be made,
        None is returned.
    :rtype: SFTPClient
    """
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
    """Upload a file to the server.

    :param sftp: Paramiko SFTPClient connection to a server.
    :type sftp: SFTPClient
    :param local_filepath: Absoluate path to file to be uploaded.
    :type local_filepath: str
    :param remote_filepath: Absoluate path to server destination.
    :type remote_filepath: str
    :returns: Indicates whether upload was successful.
    :rtype: bool
    """
    try:
        sftp.put(local_filepath, remote_filepath)
        result = True
    except:
        print("Problems encountered with the "
              "specified local or destination directories.")
        result = False
    return result
