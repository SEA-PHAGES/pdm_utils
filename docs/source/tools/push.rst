.. _push:

push: upload a database to a public server
==========================================


Data can be uploaded to the :hatfullserver:`Hatfull lab’s public server <>` using the ``pdm_utils push`` tool::

    > python3 -m pdm_utils push  -f Actinobacteriophage.sql

The '-f' flag indicates a specific file needs to be uploaded. Alternatively, a directory of files can be indicated using the '-d' flag::

    > python3 -m pdm_utils push  -d ./new_data/
