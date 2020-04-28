Running the test suites
=======================

The directory of integration and unit tests can be stored completely separately from the actual package directory. However, in order for ``unittest`` to be able to import the ``pdm_utils`` module that needs to be tested, ``unittest`` needs to be run in the directory immediately above the module. Follow the steps below to run the tests.

    1. Navigate to the src directory::

        > cd /path/to/pdm_utils/src/

    2. To run all tests::

        > python3 -m unittest discover ../tests

    3. To run only unit tests or integration tests::

        > python3 -m unittest discover ../tests/unit/
        > python3 -m unittest discover ../tests/integration/

    4. To run tests on only a specific module::

        > python3 -m unittest discover ../tests/integration/ -p test_phamerator.py


For integration tests that require a MySQL database, it is expected that MySQL user 'pdm_anon' exists with password 'pdm_anon' that has all privileges for 'pdm_test_' databases. Log in to MySQL as the 'root' user and execute the following commands::

    mysql> CREATE USER 'pdm_anon'@'localhost' IDENTIFIED BY 'pdm_anon';
    mysql> GRANT ALL PRIVILEGES ON `pdm\_test\_%`.* TO 'pdm_anon'@'localhost';
    mysql> GRANT SELECT ON *.* TO 'pdm_anon'@'localhost';
    mysql> FLUSH PRIVILEGES;


For the find_domains pipeline integration_heavy tests, it is expected that the path to the NCBI CDD is::

    ~/Databases/Cdd_LE

It may be better to store the CDD in a separate directory though, in which case a symlink can be created::

    > ln -s /path/to/CDD ~/Databases/Cdd_LE
