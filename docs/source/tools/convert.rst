.. _convert:


convert: upgrade and downgrade a database schema
================================================

Occasionally, the structure (schema) of the MySQL database is modified, which can result in certain types of data being stored differently. Even though the actual data may not have changed, these schema changes may prevent other tools from interacting with and retrieving data from the database. The ``pdm_utils convert`` tool is built to upgrade and downgrade the schema of a the database on a local computer to ensure your tools can still access the data.

To upgrade or downgrade a database (e.g. Actinobacteriophage) to the current version maintained in the Hatfull lab::

    > python3 -m pdm_utils convert Actinobacteriophage


To convert a database to a specific schema version that your tool requires, the schema version can be selected with the '-s' flag::

    > python3 -m pdm_utils convert Actinobacteriophage -s 5

Upgrading or downgrading to some schemas result in some tables or columns to be removed and added (and therefore auto-populated), resulting in some data to be lost or possibly inaccurate. After schema conversion, any tables or columns that may now have missing or inaccurate data are listed.

To avoid overwriting your primary database during conversion, a new database name can be indicated using the '-n' flag (although, if a database already exists with the same name, it will be overwritten)::

    > python3 -m pdm_utils convert Actinobacteriophage -s 5 -n Actinobacteriophage_s5

Additionally, if your local database schema doesn't exactly match one of the pre-defined schemas, MySQL may encounter a problem during conversion that results in the ``convert`` tool to terminate. This will result in a database that may only be partially converted to the intended schema. In these cases, the schema can be manually converted by executing the individual SQL statements in each conversion step. The list of MySQL conversion statements associated with each version change are stored in the :pdmutils:`pdm_utils GitHub repository <>`. Clone the repository onto your local computer. The folder of the conversion scripts is located at 'pdm_utils/schema/conversion_scripts/'.

For instance, to manually convert a schema from version 6 to 5, execute each SQL statement in 'pdm_utils/schema/conversion_scripts/downgrade_6_to_5.sql' in mysql:

    1. Log in to MySQL (enter your password when prompted)::

        > mysql -u root -p

    2. Execute each query from the conversion script, such as below::

        mysql> ALTER TABLE `gene_domain` CHANGE `QueryEnd` `query_end` int(10) unsigned NOT NULL;
