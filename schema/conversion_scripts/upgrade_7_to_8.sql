# MySQL script to upgrade the database schema from version 7 to 8.
ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_3`;
ALTER TABLE `gene`
  ADD FOREIGN KEY (`PhamID`)
  REFERENCES `pham` (`PhamID`)
  ON DELETE SET NULL;
UPDATE `version` SET `SchemaVersion` = 8;
