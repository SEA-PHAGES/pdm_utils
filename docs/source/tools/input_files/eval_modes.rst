.. _evalmodes:

Evaluation modes and evaluation flags
=====================================


Evaluation flags
----------------

Many ``import`` QC steps need to be performed on every genome (such as confirming the nucleotide sequence is not already present in the database under a separate name). However, a MySQL database may store data for diverse types of genomes, and some evaluations are dependent on factors such as the annotation status, the authorship, or the data source. As a result, some QC steps can be turned on (yes) and off (no) depending on the type of genome being imported.


check_replace
*************

Should unexpected genome replacements be reported?


import_locus_tag
****************

Should CDS feature locus_tags be imported?


check_locus_tag
***************

Should the structure of CDS feature locus_tags be evaluated?


check_description_tally
***********************

Should the number of descriptions be evaluated?


check_description_field
***********************

Should CDS descriptions in unexpected fields be reported?


check_description
*****************

Should unexpected CDS descriptions be reported?


check_trna
**********

Should tRNA features be evaluated?


check_id_typo
*************

Should genome ID typos be reported?


check_host_typo
***************

Should host typos be reported?


check_author
************

Should unexpected authors be reported?


check_gene
**********

Should the CDS 'gene' qualifier be evaluated?


check_seq
*********

Should the nucleotide sequence be evaluated?


check_coords
************

Should feature duplication be evaluated?





Evaluations modes
-----------------

In order to manage which evaluations are implemented, the evaluation mode is specified for each ticket:

draft
*****

For new, automatically-generated 'draft' annotations, since they are likely to have several types of automatically-generated errors that can be ignored.

final
*****

For new, manually reviewed 'final' annotations that are expected to have very few errors.

auto
****

For manually reviewed 'final' annotations that are automatically retrieved from GenBank. Some errors are ignored since they are already publicly available.


misc
****

For genome annotations created from an external source. Since it is not always known how they have been annotated, and since any errors may not be able to be changed, certain types of errors are ignored.

custom
******

For manual selection of the specific evaluation flags that should be performed if none of the other four preset run modes are appropriate.
