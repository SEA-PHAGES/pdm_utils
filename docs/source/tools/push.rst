.. _push:

push: upload a database to a public server
==========================================


Data can be uploaded to the :hatfullserver:`Hatfull lab’s public server <>` using the ``pdm_utils push`` tool::

    > python3 -m pdm_utils push -f Actino_Draft.sql

The '-f' flag indicates a specific file needs to be uploaded. Alternatively, a directory of files can be indicated using the '-d' flag::

    > python3 -m pdm_utils push -d ./new_data/

As with other pipelines, use of the :ref:`config_file` option can automate accessing the server to upload data.
