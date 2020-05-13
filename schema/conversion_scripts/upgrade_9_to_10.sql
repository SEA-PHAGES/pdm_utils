# MySQL script to upgrade the database schema from version 9 to 10.
ALTER TABLE `gene` MODIFY COLUMN `Translation` blob DEFAULT NULL;
UPDATE `version` SET `SchemaVersion` = 10;
