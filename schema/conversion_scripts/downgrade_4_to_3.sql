# MySQL script to downgrade the database schema from version 4 to 3.
UPDATE `gene` SET `translation` = NULL WHERE LENGTH(`translation`) > 4000;
ALTER TABLE `gene` MODIFY `translation` VARCHAR(4000) DEFAULT NULL;
DROP TABLE `tmrna`;
DROP TABLE `trna`;
DROP TABLE `trna_structures`;
UPDATE `version` SET `schema_version` = 3;

### DATA_LOSS_SUMMARY
# LOST_TABLE:tmrna
# LOST_TABLE:trna
# LOST_TABLE:trna_structures
# LOST_DATA:gene.translation
