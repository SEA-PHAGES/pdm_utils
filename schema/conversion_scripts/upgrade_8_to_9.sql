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
ALTER TABLE `trna` MODIFY COLUMN `PhageID` varchar(25) NOT NULL AFTER `TrnaID`;
ALTER TABLE `trna` MODIFY COLUMN `LocusTag` varchar(35) DEFAULT NULL AFTER `Note`;
ALTER TABLE `trna` CHANGE `TrnaID` `GeneID` varchar(35) NOT NULL;
ALTER TABLE `trna` ADD COLUMN `Name` varchar(50) NOT NULL AFTER `Length`;
ALTER TABLE `trna` DROP FOREIGN KEY `trna_ibfk_1`;
ALTER TABLE `trna`
  ADD CONSTRAINT `trna_ibfk_1`
  FOREIGN KEY (`PhageID`)
  REFERENCES `phage` (`PhageID`)
  ON DELETE CASCADE
  ON UPDATE CASCADE;

ALTER TABLE `tmrna` ADD KEY (`PhageID`);
ALTER TABLE `tmrna`
  ADD CONSTRAINT `tmrna_ibfk_1`
  FOREIGN KEY (`PhageID`)
  REFERENCES `phage` (`PhageID`)
  ON DELETE CASCADE
  ON UPDATE CASCADE;
ALTER TABLE `tmrna` MODIFY COLUMN `PhageID` varchar(25) NOT NULL AFTER `TmrnaID`;
ALTER TABLE `tmrna` MODIFY COLUMN `LocusTag` varchar(35) DEFAULT NULL AFTER `Note`;
ALTER TABLE `tmrna` CHANGE `TmrnaID` `GeneID` varchar(35) NOT NULL;
ALTER TABLE `tmrna` ADD COLUMN `Length` mediumint(9) NOT NULL AFTER `Stop`;
ALTER TABLE `tmrna` ADD COLUMN `Name` varchar(50) NOT NULL AFTER `Length`;

UPDATE `version` SET `SchemaVersion` = 9;

### DATA_LOSS_SUMMARY
# LOST_COLUMN:trna.Sequence
# LOST_COLUMN:trna.Product
# LOST_COLUMN:trna.InfernalScore
# INACCURATE_COLUMN:trna.Source
# INACCURATE_COLUMN:trna.Name
# INACCURATE_COLUMN:tmrna.Name
# INACCURATE_COLUMN:tmrna.Length
