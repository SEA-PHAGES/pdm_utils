.. _export:

export: export data from a database
===================================

This tool is used to export data from a MySQL database in a variety of formats, including:

    1. A SQL file that represents the entire database.
    2. CSV-formatted tables for selected tables.
    3. Biopython SeqIOFormatted files for selected genomes, such as:

        a. GenBank-formatted flat files (.gb)
        b. Fasta-formatted files (.fasta)

Main Pipelines
--------------

SQL File Export
_______________

Export the entire database as a SQL file::

    > python3 pdm_utils export Actinobacteriophage sql

Database information stored in a SQL database may need to be uploaded to a server for end-users.  This option allows a file to be exported from MySQL into a single file that can be easily uploaded to a server (e.g. Actinobacteriophage.sql). 
The database version is tracked as an integer in the Version field of the *version* table, and a version file is also generated (e.g. Actinobacteriophage.version), which is a text file that contains a single integer corresponding to the database version.



CSV File Export
_______________

Export specific tables from the database into CSV-formatted files::

    > python3 pdm_utils export Actinobacteriophage csv ...

Database information stored in a SQL database may need to be viewed in a more universal file.  This option exports database information from a table in a comma-separated-values file that can be viewed in common spreadsheet or text-editing softwares (e.g. phage.sql).


BioPython SeqIO Files
_____________________

Export genomes into biologically-relevant formatted files::

    > python3 pdm_utils export Actinobacteriophage gb ...

Database information stored in a SQL database may be representative of biological constructs, and common bio-file currencies are often used to exchange and update information in the SQL database.  This option exports database information to a formatted file that are common for bio-softwares (e.g. Trixie.gb, Trixie.fasta etc.)


Basic Export Options
--------------------

Changing the export folder path
_______________________________

Change the path where a new directory containing exported files will be created::

    > python3 pdm_utils export Actinobacteriophage gb -o /new/folder/path

    > python3 pdm_utils export Actinobacteriophage csv --folder_path /new/folder/path

The command-line flag **-o** or **--folder_path** followed by the path to a desired directory will set the path where a the new directory will be created.

Changing the export folder name
_______________________________

Change the path where a new directory containing exported files will be created::

    > python3 pdm_utils export Actinobacteriophage fasta -m new_folder_name

    > python3 pdm_utils export Actinobacteriophage sql --folder_name new_folder_name

The command-line flag **-m** or **--folder_name** followed by a name will set the name of the new directory to be created.

Toggling the verbosity of export 
________________________________

Toggle on export progress statements::

    > python3 pdm_utils export Actinobacteriophage csv -v

    > python3 pdm_utils export Actinobacteriophage sql --verbose 

The command-line flag **-v** or **--verbose** followed by the path to a desired directory will toggle on progress report and status statements (verbosity).

Import and Selection Export Options
-----------------------------------

Changing the table
__________________

Csv or SeqIO option to change the database table centered on for data export.::

    > python3 pdm_utils export Actinobacteriophage gb -t phage

    > python3 pdm_utils export Actinobacteriophage csv --table gene

The command-line flag **-t** or **--table** followed by a valid table from the selected MySQL database from which data is selected to be exported.  
Changing the table for csv export will change which columns are selected for export while changing the table for BioPython SeqIO file types will determine the data the formatted file will present.

Importing values with the command line
______________________________________

Csv or SeqIO option to pre-select data for export.::

    > python3 pdm_utils export Actinobacteriophage gb -in Trixie

    > python3 pdm_utils export Actinobacteriophage csv --import_names D29 L5

The command-line flag **-in** or **--import_names** followed by primary-key values from the database table selected for export (see *Changing the table*) begins export conditioned on the given set of values.

Importing values from a file
____________________________

Csv or SeqIO option to pre-select data for export.::

    > python3 pdm_utils export Actinobacteriophage gb -if /path/to/file

    > python3 pdm_utils export Actinobacteriophage csv --import_file /path/to/file

The command-line flag **-if** or **--import_file** followed by a comma-separated-values file to be read for values.  The first row of this file will be used as primary-key values from the database table selected for export (see *Changing the table*) to condition export on (similar to *Importing values with the command line*).

Including additional csv export columns
_______________________________________

Csv option to add additional columns from the database for data export.::
    
    > python3 pdm_utils export Actinobacteriophage csv -ic gene.GeneID

    > python3 pdm_utils export Actinobacteriophage csv --include_columns gene.PhamID gene.Notes

The command-line flag **-ic** or **--include_columns** followed by a MySQL-formatted column from the MySQL database selected for export to additionally be exported.
Included columns must follow the format *table*.\ *column* and can be columns from different tables than the one selected for export (see *Changing the table*).

Excluding csv export columns
____________________________

Csv option to exclude columns from the database for data export.::
    
    > python3 pdm_utils export Actinobacteriophage csv -ec phage.Subcluster

    > python3 pdm_utils export Actinobacteriophage csv --exclude_columns phage.Length

The command-line flag **-ec** or **--exclude_columns** followed by a MySQL-formatted column from the MySQL database selected for export tagged to not be exported.  
Excluded columns must follow the format *table*.\ *column*  and can be columns from different tables than the one selected for export (see *Changing the table*).

Filtering and Organization Export Options
-----------------------------------------

Filtering export
________________

Csv or SeqIO option to filter data retrieved from the database.::

    > python3 pdm_utils export Actinobacteriophage gb -f "phage.Cluster = A AND phage.Subcluster IS NOT NULL"

    > python3 pdm_utils export Actinobacteriophage csv --where "domain.Description LIKE %helix-turn-helix% OR gene.Notes = 'helix-turn-helix DNA binding protein'"

The command-line flat **-f** or **--where** followed by a MySQL-formatted WHERE expression clauses separated by ANDs and ORs.
Clauses can be expressed with the following format *table*.\ *column* *[operator]* *value* and can be using columns from different tables than the one selected for export (see *Changing the table*)

Grouping export
_______________

Csv option to exclude columns from the database for data export.::
    
    > python3 pdm_utils export Actinobacteriophage csv -g phage.Status

    > python3 pdm_utils export Actinobacteriophage csv --group_by phage.Cluster

The command-line flag **-g** or **--group_by** followed by a MySQL-formatted column from the MySQL database to group the data by for export.  Grouping creates multiple subdirectories during export, and additional groups layer the subdirectories and group within already formed groups.
Group by columns must follow the format *table*.\ *column*  and can be columns from different tables than the one selected for export (see *Changing the table*).

Sorting export
______________

Csv option to exclude columns from the database for data export.::
    
    > python3 pdm_utils export Actinobacteriophage csv -s phage.Length

    > python3 pdm_utils export Actinobacteriophage csv --order_by phage.PhageID phage.Subcluster

The command-line flag **-s** or **--order_by** followed by a MySQL-formatted column from the MySQL database to sort the data by for export.  Ordering sorts the data exported and additional orderings subsort the data.
Order by columns must follow the format *table*.\ *column* and can be columns from different tables than the one selected for export (see *Changing the table*).

Additional Export Options
-------------------------

Concatenating SeqIO files
_________________________

SeqIO option to add all export data into one contiguous formatted file.::

    > python3 pdm_utils export Actinobacteriophage gb -cc

    > python3 pdm_utils export Actinobacteriophage gb --concatenate

The command line flag **-cc** or **--concatenate** toggles the concatenation of exported SeqIO formatted flat files.

Including sequence data
_______________________

Csv option to include all sequence and translation data.::
    
    > python pdm_utils export Actinobacteriophage csv -sc

    > python pdm_utils export Actinobacteriophage csv -sequence_columns

The command line flag **-sc** or **--sequence_columns** toggles the inclusion of sequence or translation type data into the csv for export.

Conserving raw byte data
________________________

Csv option to conserve and export raw byte data.::

    > python pdm_utils export Actinobacteriophage csv -rb

    > python pdm_utils export Actinobacteriophage csv --raw_bytes

The command line flag **-rb** or **--raw_bytes** toggles off the conversion of blob and byte-type data flagged for export, exporting the raw-byte format of the data.

