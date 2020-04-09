# MySQL script to downgrade the database schema from version 9 to 8.
# Add column `Sequence` to `trna`;
ALTER TABLE `trna` ADD COLUMN `Sequence` varchar(100) NOT NULL;     # Note: Artificially restricts sequence length - may truncate some
ALTER TABLE `trna` ADD COLUMN `Product` blob;                       # Note: Will be totally empty and may not reflect the flat file
ALTER TABLE `trna` ADD COLUMN `InfernalScore` decimal(4,2);         # Note: Will be totally empty and thus not reflect reality
ALTER TABLE `trna` DROP `Source`;                                   # Note: Will lose track of which program(s) support this trna

UPDATE `trna`
  LEFT JOIN `phage` ON `trna`.`PhageID` = `phage`.`PhageID`
  SET `trna`.`Sequence` = SUBSTRING(`phage`.`Sequence` FROM `trna`.`Start` FOR `trna`.`Length`);   # Populate sequence column

CREATE TABLE `trna_structures` (
  `Sequence` varchar(100) NOT NULL,
  `Structure` varchar(300) NOT NULL,
  PRIMARY KEY (`Sequence`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

INSERT IGNORE INTO `trna_structures` SELECT `Sequence`, `Structure` FROM `trna`;    # Should copy sequences and structures to `trna_structures`
ALTER TABLE `trna` DROP `Structure`;                                                # Can now delete structure from `trna`

# Modify the set of allowable amino acids - 20 standard plus Undet, OTHER
# NOTE: Anything that was previously non-standard and non-'OTHER' will be 'OTHER' - may not be totally accurate
UPDATE `trna` SET `AminoAcid` = 'OTHER' WHERE `AminoAcid` in ('fMet', 'Ile2', 'Pyl', 'Sec', 'Stop');
ALTER TABLE `trna` MODIFY COLUMN `AminoAcid` enum('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val','Undet','OTHER') NOT NULL;