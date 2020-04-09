Schema refactoring
==================

Any changes made to the structure (schema) of the database (in the form of schema refactoring, schema improvements, and data migration) should be tracked. In order to do this, paired upgrade/downgrade scripts should be created, so that the schema changes can be implemented or reversed if needed.

    1. Determine which aspects of the schema should be changed.

    2. Create a MySQL script that contains the statements needed to make all changes (including incrementing version.SchemaVersion).

    3. In the MySQL command line utility, manually execute each statement on a test database to verify the necessary changes are successful.

    4. Once all statements are constructed, execute the entire MySQL 'upgrade' script at the command line to ensure it works properly::

        > mysql -u root -p <test database> < upgrade_script.sql

    5. Now create a 'downgrade' script to undo the changes. The order of the downgrade statements should be in the reverse order as the upgrade statements. As with the upgrade statements, first test each statement individually in the MySQL command line utility, then test the entire downgrade script at the command line.

    6. It is important that the database schema created when upgrading from an earlier schema version is identical to the database schema created when downgrading from a later schema version. Using a test database, this can be determined with the new paired upgrade/downgrade scripts as follows:

        1. Export the empty schema before upgrade::

            > mysqldump --no-data -u root -p --skip-comments <db_name> > db_schema_before.sql

        2. Run the upgrade script::

            > mysql -u root -p <test database> < upgrade_script.sql

        3. Run the downgrade script::

            > mysql -u root -p <test database> < downgrade_script.sql

        4. Export the empty schema after downgrade::

            > mysqldump --no-data -u root -p --skip-comments <db_name> > db_schema_after.sql

        5. Check the difference between the empty schemas. Other than AUTO-INCREMENT values, there should be no substantial differences::

            > diff db_schema_before.sql db_schema_after.sql

        6. If the conversion round-trip does not produce an identical empty schema, modify the upgrade or downgrade statements accordingly.

    7. Incorporate the upgrade and downgrade statements into the ``pdm_utils`` schema_conversions module so that they can be implemented using the Python package.

    8. In the convert module, edit the CURRENT_VERSION and MAX_VERSION variables accordingly.

    9. Use the convert module to upgrade the primary production database to the new schema version. This will convert the schema and update version.SchemaVersion.

    10. A history of each unique database schema is stored under /misc/schemas/. Create an empty schema of the upgraded database::

        > mysqldump --no-data -u root -p --skip-comments <db_name> > db_schema_<new schema int>.sql

    11. Add the sql file to the schemas directory in the repo.

    12. Update the schema_updates.txt history file with the changes, including a summary of the types of changes implemented.

    13. Generate a schema map and update all sections of the user guide (see below).

        1. Open MySQL Workbench and connect to the server.

        2. Under the Database menu, select Reverse Engineer.

        3. Choose the database of interest.

        4. Manually move table icons so they are intuitively arranged.

        5. Under File, select Export, then select Export as Single Page PDF.

        6. Open the PDF in Preview, and under File, select Export, then select Format JPEG 300 resolution.

        7. Add the JPEG to the repo in the user guide directory.

    14. Update the user guide as needed with information about the new schema:

        - page describing the current database
        - page describing prior schema version schema maps
        - page describing schema version changelog

    15. Update files requires for integration tests:

        1. Copy the new empty schema file to the following directory::

            tests/integration/test_files/

        2. Rename the file to 'test_schema_<schema_version>.sql'
        3. Remove the test schema file of the prior version.
