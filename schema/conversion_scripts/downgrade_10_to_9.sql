# MySQL script to downgrade the database schema from version 10 to 9.
UPDATE `gene` SET `Translation` = NULL WHERE LENGTH(`Translation`) > 5000;
ALTER TABLE `gene` MODIFY COLUMN `Translation` varchar(5000) DEFAULT NULL;
UPDATE `version` SET `SchemaVersion` = 9;

### DATA_LOSS_SUMMARY
# LOST_DATA:gene.Translation
