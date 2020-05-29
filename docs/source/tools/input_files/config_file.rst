.. _config_file:

Database management configuration file
======================================

Every pipeline can accept a configuration file that contains user-specific information to enable the pipeline to be run automatically, including MySQL login credentials, NCBI credentials (for rapid retrieval of GenBank records), and server information and credentials (to retrieve data from, or upload data to, a server).

The configuration file is formatted to meet requirements for the :configparser:`configparser <>` module. Below is an example of how it is structured with a brief description of each parameter::

    [ncbi]                      # Settings for connecting to NCBI server
    email=jane@acme.com         # Your email address
    api_key=123456789           # Your NCBI api key
    tool=DataRetrievalTool      # The name of your script

    [mysql]                     # Settings for connecting to local MySQL
    user=janedoe                # Your MySQL login username
    password=random             # Your MySQL login password

    [upload_server]             # Settings to upload file(s) to a server
    host=abc.def.com            # Name of the host server to upload file(s)
    dest=/path/to/folder/       # Server path to where the file(s) will be uploaded
    user=janedoe                # Your server login username
    password=random             # Your server login password

    [download_server]           # Settings to download file(s) from a server
    url=http://abc.def.com/     # URL hosting the file to download


Not every pipeline requires all information in the config file, so in general, individual sections or rows are optional. If a pipeline needs information not present in the config file, it will prompt the user for the requisite information. For NCBI, if you do not provide an email address or api key, your GenBank queries may be throttled or blocked (please refer to :ncbi:`NCBI <>` for more details). An example config file is available on the ``pdm_utils`` :pdmutils:`source code repository <>`.
