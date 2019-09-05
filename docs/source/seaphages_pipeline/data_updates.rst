Data updates
============


Although it is recommended that genome and gene data are imported and updated directly from flat files using the import script, sometimes it may be necessary to modify or update specific fields for specific phages. This can be accomplished using a collection of update_[field].py scripts, in which the [field] refers to the specific field that needs to be updated [Table 8].

Each script is similarly structured, is pointed to the appropriate table in the database to implement updates, and requires two arguments [Table 9]. In contrast to the 11-field import tickets, a list of simple update tickets are stored in a csv-formatted update_table containing two columns. The first column contains the primary key in the database table for which the target field will be updated, and the second column contains the new data that will populate the target field. For example, in the update_gene_description.py script, gene descriptions for hundreds of genes from hundreds of phages can be updated at once; the update_table contains a list of unique GeneIDs that are derived from the Gene table and accompanying gene descriptions, and the script implements these changes in the Gene table. Unlike import tables though, the csv-formatted update tables need to be manually constructed.

It is also important to note that this collection of scripts is under-developed. They do not perform any QC checks on the data in the file. For instance, if a GeneID from the update_table is not found in the Gene table in the database, the ticket fails to be implemented but no warning or error is issued. Therefore, although these scripts are valuable, it is important to use caution when utilizing them.
