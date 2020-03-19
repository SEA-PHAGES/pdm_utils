How to create the database for integration tests
================================================

Some integration tests require a pre-defined, non-empty test database, called 'test_db.sql'. This database is created from Actinobacteriophage by deleting all but 15 phages:


1.  Drop the test_db database (if it exists)::

        mysql -u root -p --execute "DROP DATABASE test_db;"

2.  Create an empty test_db database::

        mysql -u root -p --execute "CREATE DATABASE test_db;"


3.  Download the most recent version of Actinobacteriophage.sql::

        python3 -m pdm_utils get_db Actinobacteriophage ./ -d

4.  Source the new Actinobacteriophage.sql into test_db::

        mysql -u root -p test_db < Actinobacteriophage.sql

5.  Log into mysql and run the following commands::

        USE test_db;
        DELETE FROM phage WHERE PhageID NOT IN ('Alice', 'Atlantean', 'Myrna', 'MichelleMyBell', 'BaxterFox', 'Octobien14', 'Aubergine', 'Lucky3', 'Constance', 'Mufasa8', 'Yvonnetastic', 'Et2Brutus', 'D29', 'L5', 'Trixie', 'Sparky');
        DELETE FROM pham WHERE PhamID NOT IN (SELECT DISTINCT PhamID FROM gene);
        DELETE FROM domain WHERE HitID NOT IN (SELECT DISTINCT HitID FROM gene_domain);
        COMMIT;
        EXIT;

6.  Dump the new database::

        mysqldump -u root -p --skip-comments test_db > test_db.sql

7.  Copy test_db.sql to the test directory::

        cp ./test_db.sql <path to repo>/pdm_utils/tests/test_files


.. csv-table:: Phages in the test database
    :file: test_db_phages.csv
