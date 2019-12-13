Misc notes that need to be added somewhere:


Should probably add in section that highlights the weekly routine database updates.
The import script is designed to handle diverse types of tickets present in a single import_table. However, the retrieve_data.py script creates separate staged directories and import tables for different types of data to be imported to minimize potential ticket conflicts. When the import_script.py is executed following the retrieve_data.py script, it is recommended that the script is executed separately for each ticket type, and in the following order: metadata updates, auto-annotated genomes from PECAAN, new preliminary final annotations from PhagesDB, auto-updated SEA-PHAGES final annotations from GenBank, and other miscellaneous tickets that need to be implemented.
