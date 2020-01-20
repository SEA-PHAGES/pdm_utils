# MySQL script to downgrade the database schema from version 8 to 7.
ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_3`;
ALTER TABLE `gene`
  ADD FOREIGN KEY (`PhamID`)
  REFERENCES `pham` (`PhamID`);
UPDATE `version` SET `SchemaVersion` = 7;
