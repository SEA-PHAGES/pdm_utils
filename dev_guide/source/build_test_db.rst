How to create the database for integration tests
================================================

Some integration tests require a pre-defined, non-empty test database, called 'test_db.sql'. This database is created from Actinobacteriophage by deleting all but 15 phages:

1.  From the command line, drop test_db and re-create it::

        mysql -u root -p --execute "DROP DATABASE test_db; CREATE DATABASE test_db;"

2.  From the command line, source Actinobacteriophage.sql into test_db::

        mysql -u root -p test_db < Actinobacteriophage.sql

3.  Log into mysql and run the following commands::

        USE test_db;
        DELETE FROM phage WHERE PhageID NOT IN ('Alice', 'Atlantean', 'Myrna', 'MichelleMyBell', 'BaxterFox', 'Octobien14', 'Aubergine', 'Lucky3', 'Constance', 'Mufasa8', 'Yvonnetastic', 'Et2Brutus', 'D29', 'L5', 'Sparky');
        COMMIT;
        EXIT;

4.  Dump the new database::

        mysqldump -u root -p --skip-comments test_db > test_db.sql

5.  Copy test_db.sql to the test directory::

        cp ./test_db.sql /pdm_utils/tests/test_files


.. csv-table:: Phages in the test database
    :file: test_db_phages.csv
