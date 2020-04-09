# MySQL script to upgrade the database schema from version 8 to 9.
ALTER TABLE `trna` ADD COLUMN `Structure` varchar(300) NOT NULL;
UPDATE `trna`
  LEFT JOIN `trna_structures` ON `trna`.`Sequence` = `trna_structures`.`Sequence`
  SET `trna`.`Structure` = `trna_structures`.`Structure`;
ALTER TABLE `trna` MODIFY COLUMN `Structure` blob DEFAULT NULL;
DROP TABLE `trna_structures`;
UPDATE `trna` SET `AminoAcid` = 'OTHER' WHERE `AminoAcid` = 'Undet';
ALTER TABLE `trna` MODIFY COLUMN `AminoAcid` enum('Ala','Arg','Asn','Asp','Cys','fMet','Gln','Glu','Gly','His','Ile','Ile2','Leu','Lys','Met','Phe','Pro','Pyl','SeC','Ser','Thr','Trp','Tyr','Val','Stop','OTHER') NOT NULL;
ALTER TABLE `trna` DROP COLUMN `Sequence`;
ALTER TABLE `trna` DROP COLUMN `Product`;
ALTER TABLE `trna` DROP COLUMN `InfernalScore`;
ALTER TABLE `trna` ADD COLUMN `Source` enum('aragorn', 'trnascan', 'both') DEFAULT NULL;


### DATA_LOSS_SUMMARY
# LOST_TABLE:trna_structures
# LOST_COLUMN:trna.Sequence
# LOST_COLUMN:trna.Product
# LOST_COLUMN:trna.InfernalScore
# INACCURATE_COLUMN:trna.Source
