# MySQL script to downgrade the database schema from version 9 to 8.
ALTER TABLE `trna` DROP COLUMN `Source`;
ALTER TABLE `trna` ADD COLUMN `InfernalScore` decimal(4,2) DEFAULT NULL AFTER `Anticodon`;
ALTER TABLE `trna` ADD COLUMN `Product` blob AFTER `Orientation`;
ALTER TABLE `trna` ADD COLUMN `Sequence` varchar(100) DEFAULT NULL AFTER `Orientation`;
UPDATE `trna`
  LEFT JOIN `phage` ON `trna`.`PhageID` = `phage`.`PhageID`
  SET `trna`.`Sequence` = SUBSTRING(`phage`.`Sequence` FROM `trna`.`Start` FOR `trna`.`Length`);
ALTER TABLE `trna` MODIFY COLUMN `Sequence` varchar(100) NOT NULL;
UPDATE `trna` SET `AminoAcid` = 'OTHER' WHERE `AminoAcid` in ('fMet', 'Ile2', 'Pyl', 'Sec', 'Stop');
ALTER TABLE `trna` MODIFY COLUMN `AminoAcid` enum('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val','Undet','OTHER') NOT NULL;

CREATE TABLE `trna_structures` (
  `Sequence` varchar(100) NOT NULL,
  `Structure` varchar(300) NOT NULL,
  PRIMARY KEY (`Sequence`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
ALTER TABLE `trna` MODIFY COLUMN `Structure` varchar(300) NOT NULL;
INSERT IGNORE INTO `trna_structures` SELECT `Sequence`, `Structure` FROM `trna`;
ALTER TABLE `trna` DROP COLUMN `Structure`;

### DATA_LOSS_SUMMARY
# INACCURATE_COLUMN:trna.Product
# INACCURATE_COLUMN:trna_structures.Sequence
# INACCURATE_COLUMN:trna.InfernalScore
# INACCURATE_COLUMN:trna.AminoAcid
# LOST_COLUMN:trna.Source
