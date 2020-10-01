.. _revise:

revise: revise data inconsistencies
===================================

This tool is to be used in combination with the ``pham_review`` pipeline to generate files formatted to suit a GenBank automated gene product annotation resubmission pipeline.

Specifically, the local ``revise`` pipeline reads in edits made to a review csv spreadsheet and searches through an indicated database to find instances of genes selected for review with product annotations dissimilar to the indicated product annotation.  The ``revise`` pipeline marks these, and generates a formatted file of the changes required to correct them.

To export the resubmission csv file from a edited review file::

    > python3 -m pdm_utils revise Actinobacteriophage local <path/to/review/file>

To run the ``revise`` pipeline in production mode, which takes into account the various flags used to denote submitted GenBank files, utilize the production flag at the command line::

    >python3 -m pdm_utils revise Actinobacteriophage local <path/to/review/file> --production

Like many other pipelines that involve exporting data in the ``pdm_utils`` package, the range of the database entries inspected by the ``revise`` pipeline can be modified with command line ``filter`` module implementations.  Filtering, grouping, and sorting can be done in the same manner as described in the ``export`` pipeline documentation::

    > python3 pdm_utils revise Actinobacteriophage local FunctionReport.csv -w "gene.Cluster NOT IN ('A', 'B', 'K')" -g phage.Cluster -s phage.PhageID

The local ``revise`` pipeline can also translate review formatted data into update ticket tables that can be used to update the database::

    >python3 pdm_utils revise Actinobacteriophage local FunctionReport.csv -ft ticket

The remote ``revise`` pipeline retrieves data from GenBank in five-column feature table format and searches through an indicated database to find discrepancies between the product annotations and starts in the local database and those stored at GenBank.  The ``revise`` pipeline marks these, and edits the retrieved GenBank files to generate five-column feature tables with feature data consistent with the local data.

To generate the revised five-column feature table format files::

    > python3 -m pdm_utils revise Actinobacteriophage remote 

And again, the remote ``revise`` pipeline can be modified with command line ``filter`` module implementations in the same manner as described in the ``export`` pipeline documentation::

    > python3 pdm_utils revise Actinobacteriophage remove -w "gene.Subcluster IN ('K1', 'K2', 'K6')
