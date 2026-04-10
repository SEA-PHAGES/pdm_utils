# MySQL script to upgrade the database schema from version 10 to 11.
ALTER TABLE `gene` MODIFY COLUMN `MembraneStatus` tinyint(1) DEFAULT '1' NOT NULL;
UPDATE `version` SET `SchemaVersion` = 10;