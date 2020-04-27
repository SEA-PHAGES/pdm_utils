# MySQL script to downgrade the database schema from version 9 to 8.
ALTER TABLE `tmrna` DROP COLUMN `Name`;
ALTER TABLE `tmrna` DROP COLUMN `Length`;
ALTER TABLE `tmrna` CHANGE `GeneID` `TmrnaID` varchar(35) NOT NULL;
ALTER TABLE `tmrna` MODIFY COLUMN `LocusTag` varchar(35) DEFAULT NULL AFTER `PhageID`;
ALTER TABLE `tmrna` MODIFY COLUMN `TmrnaID` varchar(35) NOT NULL AFTER `PhageID`;
ALTER TABLE `tmrna` DROP FOREIGN KEY `tmrna_ibfk_1`;
ALTER TABLE `tmrna` DROP KEY `PhageID`;

ALTER TABLE `trna` DROP FOREIGN KEY `trna_ibfk_1`;
ALTER TABLE `trna`
  ADD CONSTRAINT `trna_ibfk_1`
  FOREIGN KEY (`PhageID`)
  REFERENCES `phage` (`PhageID`);
ALTER TABLE `trna` DROP COLUMN `Name`;
ALTER TABLE `trna` CHANGE `GeneID` `TrnaID` varchar(35) NOT NULL;
ALTER TABLE `trna` MODIFY COLUMN `LocusTag` varchar(35) DEFAULT NULL AFTER `PhageID`;
ALTER TABLE `trna` MODIFY COLUMN `TrnaID` varchar(35) NOT NULL AFTER `PhageID`;
ALTER TABLE `trna` DROP COLUMN `Source`;
ALTER TABLE `trna` ADD COLUMN `InfernalScore` decimal(4,2) DEFAULT NULL AFTER `Anticodon`;
ALTER TABLE `trna` ADD COLUMN `Product` blob AFTER `Orientation`;
ALTER TABLE `trna` ADD COLUMN `Sequence` varchar(100) DEFAULT NULL AFTER `Orientation`;
UPDATE `trna`
  LEFT JOIN `phage` ON `trna`.`PhageID` = `phage`.`PhageID`
  SET `trna`.`Sequence` = SUBSTRING(`phage`.`Sequence` FROM (`trna`.`Start` + 1) FOR `trna`.`Length`);
ALTER TABLE `trna` MODIFY COLUMN `Sequence` varchar(100) NOT NULL;
UPDATE `trna` SET `AminoAcid` = 'OTHER' WHERE `AminoAcid` in ('fMet', 'Ile2', 'Pyl', 'Sec', 'Stop');
ALTER TABLE `trna` MODIFY COLUMN `AminoAcid` enum('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val','Undet','OTHER') NOT NULL;
CREATE TABLE `trna_structures` (
  `Sequence` varchar(100) NOT NULL,
  `Structure` varchar(300) NOT NULL,
  PRIMARY KEY (`Sequence`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
UPDATE `trna` SET `Structure` = '' WHERE `Structure` is NULL;
UPDATE `trna` SET `Structure` = SUBSTRING(`Structure` FROM 1 FOR 300);
ALTER TABLE `trna` MODIFY COLUMN `Structure` varchar(300) NOT NULL;
INSERT IGNORE INTO `trna_structures` SELECT `Sequence`, `Structure` FROM `trna`;
ALTER TABLE `trna` DROP COLUMN `Structure`;

UPDATE `version` SET `SchemaVersion` = 8;

### DATA_LOSS_SUMMARY
# INACCURATE_COLUMN:trna.Product
# INACCURATE_COLUMN:trna_structures.Sequence
# INACCURATE_COLUMN:trna.InfernalScore
# INACCURATE_COLUMN:trna.AminoAcid
# LOST_COLUMN:trna.Source
# LOST_COLUMN:trna.Name
# LOST_COLUMN:tmrna.Name
# LOST_COLUMN:tmrna.Length
