How to create the database for integration tests
================================================

Some integration tests require a pre-defined, non-empty test database, called 'test_db_filled.sql'. This database is created from Actino_Draft by deleting all but 16 phages:


1.  Drop the pdm_test_db database (if it exists)::

        mysql -u root -p --execute "DROP DATABASE pdm_test_db;"

2.  Create an empty pdm_test_db database::

        mysql -u root -p --execute "CREATE DATABASE pdm_test_db;"


3.  Download the most recent version of Actino_Draft.sql::

        python3 -m pdm_utils get_db Actino_Draft server -d -o ./

4.  Source the new Actino_Draft.sql into pdm_test_db::

        mysql -u root -p pdm_test_db < Actino_Draft.sql

5. Convert the database to the correct schema, if needed::

        python3 -m pdm_utils convert pdm_test_db

6.  Log into mysql and run the following commands::

        USE pdm_test_db;
        DELETE FROM phage WHERE PhageID NOT IN ('Alice', 'Atlantean', 'Myrna', 'MichelleMyBell', 'BaxterFox', 'Octobien14', 'Aubergine', 'Lucky3', 'Constance', 'Mufasa8', 'Yvonnetastic', 'Et2Brutus', 'D29', 'L5', 'Trixie', 'Sparky');
        DELETE FROM pham WHERE PhamID NOT IN (SELECT DISTINCT PhamID FROM gene);
        DELETE FROM domain WHERE HitID NOT IN (SELECT DISTINCT HitID FROM gene_domain);
        COMMIT;
        EXIT;

7.  Dump the new database::

        mysqldump -u root -p --skip-comments pdm_test_db > test_db_filled.sql

8.  Copy test_db_filled.sql to the test directory::

        cp ./test_db_filled.sql <path to repo>/pdm_utils/tests/test_files


.. csv-table:: Phages in the test database
    :file: test_db_phages.csv
