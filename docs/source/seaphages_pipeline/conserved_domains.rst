Identifying conserved domains
=============================


The NCBI maintains a conserved domain database (CDD) in which the functions of protein sequences are systematically categorized and organized [ref Marchler-Bauer 2011]. Every gene product in PhameratorDB is evaluated for conserved domains using a local copy of the CDD [Figure 1], and this conserved domain data is stored in the database using the cdd_pp.py script. The script requires two arguments [Table 10].

In the Gene table, there is a field called Cdd_Status. When new phage genomes are added to PhameratorDB, the Cdd_Status field for each new gene is set to 0. The cdd_script.py script retrieves gene products (stored in the Translation field of the Gene table) for all genes with cdd_status < 1. The rpsblast+ package is used to identity conserved domains using BLAST with an e-value threshold = 0.001. For each gene, retrieved CDD data is inserted into the Domain and Gene_Domain tables, and the Cdd_Status field in the Gene table is set to 1 so that this gene is not re-processed during subsequent rounds of updates.
