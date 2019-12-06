# MySQL script to downgrade Phamerator database schema from version 4 to 3.
# Note: Data in several tables will be lost.
# Note: some data in gene.translation will be lost.
UPDATE `gene` SET `translation` = NULL WHERE LENGTH(`translation`) > 4000;
ALTER TABLE `gene` MODIFY `translation` VARCHAR(4000) DEFAULT NULL;
DROP TABLE `tmrna`;
DROP TABLE `trna`;
DROP TABLE `trna_structures`;
UPDATE `version` SET `schema_version` = 3;
