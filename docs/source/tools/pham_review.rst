.. _review:

review: review data for consistency
===================================

This tool is used to generate files used to review and revise gene product annotations, focusing on genes with similar amino acid sequences that have discrepant product annotations.

Specifically, the review pipeline focuses on genes within the same pham that have discrepant product annotations in order to review and revise product annotation inaccuracies that are due to human error or outdated data.

To export a csv spreadsheet of phams with discrepant product annotations and other relevant data from a database::

    > python3 -m pdm_utils pham_review Actinobacteriophage

Like many other pipelines that involve exporting data in the ``pdm_utils`` package, the output of these phams can be modified with command line ``filter`` module implementations.  Filtering, grouping, and sorting can be done in the same manner as described in the ``export`` module docs::

    > python3 -m pdm_utils pham_review Actinobacteriophage -w "phage.Cluster = 'A'" -g phage.Subcluster -s gene.Name

For attempting to revise specific phams that do not necessarily have genes with discrepant product annotations, the process of 'review' can be toggled off with the command-line flag **-r** or **--review**::

    > python3 -m pdm_utils pham_review Actinobacteriophage -nr


References to the database that the review was performed on, and the profile of the phages of the genes within the phams selected for review can be obtained by generating a summary report with the command-line flag **-sr** or **--summary_report**::

    > python3 -m pdm_utils pham_review Actinobacteriophage -sr

In addition, other supplementary information files that might help with the review process can be generated.  Translation data and data about the respective phage genomes are available in the gene reports generated with the command-line flag **-gr** or **--gene_reports**::

    > python3 -m pdm_utils pham_review Actinobacteriophage -gr

A comprehensive profile of the pham including Conserved Domain Database data and the most common annotations of adjacent genes can be generated using the command-line flag **-psr** or **--pham_summary_reports**::

    > python3 -m pdm_utils pham_review Actinobacteriophage -psr

A complete review with all reports included can be done with the command-line flag **-a** or **--all_reports**::

    > python3 -m pdm_utils pham_review Actinobacteriophage -a

The ``pham_review`` pipeline can also be targetted for gene data in the Actinobacteriophage database that has been recorded as submitted to GenBank and accessible for change as part of a pipeline of exchanging/updating data with GenBank::
    
    > python3 -m pdm_utils pham_review Actinobacteriophage -p
