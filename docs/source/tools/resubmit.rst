.. _resubmit:

resubmit: submit data updates to GenBank
========================================

This tool is to be used in combination with the review pipeline to generate files 
formatted to suit a GenBank automated gene product annotation resubmission pipeline.  

Specifically, the resubmit pipeline reads in edits made to a review csv spreadsheet 
and searches through an indicated database to find instances of genes selected for
 review with product annotations dissimilar to the indicated product annotation.  The
 resubmit pipeline marks these, and generates a formatted file of the changes required 
 to correct them.

To export the resubmission csv file from a edited review file::

    > python3 -m pdm_utils resubmit Actinobacteriophage <path/to/review/file>

Like many other pipelines that involve exporting data in the ``pdm_utils`` package, 
the range of the database entries inspected by the resubmit pipeline can be modified by
 with command line ``pdm_utils filter`` module implementations.  Filtering, grouping, 
 and sorting can be done in the same manner as described in the ``pdm_utils export_db``
 module docs::
    
    > python3 pdm_utils resubmit Actinobacteriophage FunctionReport.csv -f "gene.PhamID='10000'" -g phage.Cluster -s phage.PhageID


