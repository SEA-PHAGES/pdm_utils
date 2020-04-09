# MySQL script to upgrade the database schema from version 8 to 9.
ALTER TABLE `trna` ADD COLUMN `Structure` blob DEFAULT NULL;

# Copy from trna_structures.Structure to trna.Structure, then delete now-defunct trna_structures
UPDATE `trna`
  LEFT JOIN `trna_structures` ON `trna`.`Sequence` = `trna_structures`.`Sequence`
  SET `trna`.`Structure` = `trna_structures`.`Structure`;
DROP TABLE `trna_structures`;

# Modify the set of allowable amino acids - 20 standard plus SeC, Pyl, generic Stop, fMet, Ile2, OTHER
# NOTE: Anything that was previously 'Undet' will be 'OTHER'; anything that was 'OTHER' may not be as specific as would now be allowed
UPDATE `trna` SET `AminoAcid` = 'OTHER' WHERE `AminoAcid` = 'Undet';
ALTER TABLE `trna` MODIFY COLUMN `AminoAcid` enum('Ala','Arg','Asn','Asp','Cys','fMet','Gln','Glu','Gly','His','Ile','Ile2','Leu','Lys','Met','Phe','Pro','Pyl','SeC','Ser','Thr','Trp','Tyr','Val','Stop','OTHER') NOT NULL;

# Delete sequence, product, and Infernal columns because they are unneeded, redundant, and non-universal (respectively)
ALTER TABLE `trna` DROP `Sequence`;
ALTER TABLE `trna` DROP `Product`;
ALTER TABLE `trna` DROP `InfernalScore`;

# Add a column to store which program(s) predict the annotated tRNA
# NOTE: This column will have all NULL values after upgrade
ALTER TABLE `trna` ADD COLUMN `Source` enum('aragorn', 'trnascan', 'both') DEFAULT NULL;