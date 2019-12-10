.. _ncbicreds:

NCBI credentials file
=====================

This file can store information about:

    1. Your email address
    2. An NCBI api key
    3. The name of your script

Each row is optional, but if you do not provide an email address or api key, your GenBank queries may be throttled or blocked (please refer to :ncbi:`NCBI <>` for more details). Each row of the credentials file needs to be specifically structured as follows::

    ncbi_email=<your.address@.email.com>
    ncbi_api_key=<36 digit key>
    ncbi_tool=<name of your pipeline>
