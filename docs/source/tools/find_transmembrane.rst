.. _find_transmembrane:

find_transmembrane: find transmembrane domains
=========================================


Several tools exist to detect domains associated with the cell membrane, including SOSUI and DeepTMHMM.  Every gene product in the database can be evaluated for the presence of signal or transmembrane domains by running one of the aforementioned software.  After running this software, the domain data is stored in the database using the ``pdm_utils find_transmembrane`` tool.

To identify membrane domains::

    > python3 -m pdm_utils find_transmembrane Actino_Draft

The argument 'Actino_Draft' indicates the name of the database that will be evaluated.

In the *gene* table, there is a field called MembraneStatus.  When new phage genomes are added, the MembraneStatus field for each new gene is set to '0'.  The ``find_transmembrane`` tool retrieves gene products (stored in the Translation field of the *gene* table) for all genes with MembraneStatus M '1'. Either SOSUI or DeepTMHMM is used to identify membrane-associated domains, and the MembraneStatus field in the *gene* table is set to 1 so that this gene is not re-processed during subsequent rounds of updates.

As with other pipelines, use of the :ref: `config_file` option can automate accessing MySQL.

By default, the ``find_transmembrane`` pipeline uses DeepTMHMM and utilizes the resources available with `Biolib <https://dtu.biolib.com/DeepTMHMM>`_. To use local resources, you can use the --run_machine flag::

    > python3 -m pdm_utils find_transmembrane Actino_Draft --run_machine local

To reset the transmembrane domain information in a database, specify the --reset flag::

    > python3 -m pdm_utils find_transmembrane Actino_Draft --reset
