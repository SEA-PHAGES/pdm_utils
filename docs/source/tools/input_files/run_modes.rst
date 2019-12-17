.. _runmodes:

Run modes and evaluation flags
==============================



Optional QC steps
-----------------

Many QC steps in the import script need to be performed on every genome (such as confirming the nucleotide sequence is not already present in the database under a separate name). However, a MySQL database may store data for diverse types of genomes, and some QC steps are dependent on factors such as the annotation status, the authorship, or the data source. As a result, some QC steps can be turned on (yes) and off (no) depending on the type of genome being imported.


use_basename
************

By default, phage names in flat files are expected to be in the SOURCE-ORGANISM field. When this QC option is selected, the name of the file (without the file extension) is used as the phage name (yes = filename is used). This option is useful when importing non-SEA-PHAGES genomes.

custom_gene_id
**************

By default, the GeneID is derived from the LOCUS_TAG. When this QC option is selected, the GeneID is created by concatenating the PhageID and CDS count (yes = the GeneID is created by concatenation). This option is useful when importing non-SEA-PHAGES genomes.

ignore_gene_id_typo
*******************

By default, a warning is issued if a GeneID does not contain the phage name, indicating there is likely a typo in the GeneID. When this QC option is selected, this warning is silenced (yes = GeneID spelling is ignored). This option is useful when importing genomes from GenBank; since the GenBank LOCUS_TAG cannot be changed, there is no need for the script to issue warnings.

ignore_description_field_check
******************************

By default, a warning is issued if gene descriptions appear to be present in fields other than the field indicated by the import ticket. When this QC option is selected, this warning is silenced (yes = import gene description data from the indicated ticket field without checking other fields). This option is useful when auto-updating SEA-PHAGES genomes from GenBank, which have been systematically annotated with descriptions in the PRODUCT field.

ignore_replace_warning
**********************

By default, a warning is issued if a genome with final status is about to be replaced with a new genome. When this QC option is selected, this warning is silenced (yes = final status is ignored). This option is useful when importing genomes from GenBank, when it is expected that final status genomes will be replaced.

ignore_trna_check
*****************

By default, tRNA features are evaluated for quality, and warnings are issued when problems are encountered. When this QC option is selected, these warnings are silenced (yes = tRNA QC is ignored). This option is useful when importing draft status genomes or genomes from GenBank.

ignore_locus_tag_import
***********************

By default, data from the GenBank LOCUS_TAG field is stored in the Gene table LocusTag field. However, the LocusTag field should only reflect data from official GenBank records. When this option is selected, LOCUS_TAG data is not imported (yes = locus tags are ignored). This option is useful when importing any genome that has not been obtained from GenBank.

ignore_phage_name_typos
***********************

By default, a warning is issued if any of the various phage name fields in the flat file contain phage name typos. When this option is selected, the warning is silenced (yes = phage name typos are ignored). This option is useful when importing non-SEA-PHAGES genomes from GenBank.

ignore_host_typos
*****************

By default, a warning is issued if any of the various host name fields in the flat file contain host name typos. When this option is selected, the warning is silenced (yes = host genus typos are ignored). This option is useful when importing non-SEA-PHAGES genomes from GenBank.

ignore_generic_author
*********************

By default, a warning is issued if the author field in the flat file contains a generic author “Lastname, Firstname”, which can be inadvertently added during genome annotation. When this option is selected, the warning is silenced (yes = generic authors are ignored). This option is useful when importing draft status genomes, or genomes from GenBank.

ignore_description_check
************************

By default, a warning is issued if gene descriptions appear to contain errors (although, this QC step is currently under-developed). When this option is selected, the warning is silenced (yes = gene description errors are ignored). This option is useful when importing draft status genomes or genomes from GenBank.



Run modes
---------

In order to manage which optional QC steps are implemented, run modes have been created that are specified for each ticket:

    - pecaan: automatically-generated draft annotations.
    - phagesdb: SEA-PHAGES final annotations retrieved from PhagesDB.
    - ncbi_auto: SEA-PHAGES final annotations retrieved from GenBank.
    - ncbi_misc: non-SEA-PHAGES annotations retrieved from GenBank.
    - custom: enables the database administrator to manually select which of the 11 QC steps should be performed if none of the other four preset run modes are appropriate.
