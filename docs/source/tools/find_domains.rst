.. _findcdd:


find_domains: find NCBI conserved domains
=========================================


The NCBI maintains a :cdd:`conserved domain database (CDD) <>` in which the functions of protein sequences are systematically categorized and organized [ref Marchler-Bauer 2011]. Every gene product in the database can be evaluated for conserved domains using a local copy of the CDD, and this conserved domain data is stored in the database using the ``pdm_utils find_domains`` tool.

To identify conserved domains::

    > python3 -m pdm_utils find_domains Actinobacteriophage /path/to/CDD/

The argument 'Actinobacteriophage' indicates the name of the database that will be evaluated. The '/path/to/CDD/' indicates where the NCBI conserved domain database is stored on your local computer.

In the *gene* table, there is a field called DomainStatus. When new phage genomes are added, the DomainStatus field for each new gene is set to '0'. The ``find_domains`` tool retrieves gene products (stored in the Translation field of the *gene* table) for all genes with DomainStatus < '1'. As part of the :blastplus:`BLAST+ package <>`, the rpsblast+ tool is used to identity conserved domains using BLAST with an e-value threshold = 0.001. For each gene, retrieved CDD data is inserted into the *domain* and *gene_domain* tables, and the DomainStatus field in the *gene* table is set to 1 so that this gene is not re-processed during subsequent rounds of updates.

``find_domains`` tries to insert all hits from rpsblast+ that pass the threshold.
However, rpsblast+ may report multiple hits to a domain within the same translation that pass the threshold. Duplicate hits are not permitted in the database, so when there is an attempt to insert a duplicated hit into the database, the user is notified::

    Warning: (1062, "Duplicate entry 'gnl|CDD|334100' for key 'hit_id'")

This message is simply a warning, and no action needs to be taken.

As with other pipelines, use of the :ref:`config_file` option can automate accessing MySQL.
