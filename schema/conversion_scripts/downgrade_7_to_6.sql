# MySQL script to downgrade the database schema from version 7 to 6.
RENAME TABLE `pham` TO `pham_color`;
CREATE TABLE `pham` (
  `GeneID` varchar(35) NOT NULL DEFAULT '',
  `Name` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`GeneID`),
  KEY `name_index` (`Name`),
  CONSTRAINT `pham_ibfk_1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_3`;
INSERT INTO `pham` (`GeneID`, `Name`)
    SELECT `GeneID`,`PhamID` FROM `gene`;
ALTER TABLE `gene` DROP `PhamID`;
ALTER TABLE `gene` DROP `Parts`;
ALTER TABLE `gene` ADD COLUMN `ID` int(10) unsigned NOT NULL AFTER `Notes`;
ALTER TABLE `gene` ADD KEY `id` (`ID`);
ALTER TABLE `gene` MODIFY `ID` int(10) unsigned NOT NULL AUTO_INCREMENT;
ALTER TABLE `pham_color` DROP PRIMARY KEY;
ALTER TABLE `pham_color` ADD COLUMN `ID` int(10) unsigned NOT NULL FIRST;
ALTER TABLE `pham_color` ADD KEY (`ID`);
ALTER TABLE `pham_color` MODIFY `ID` int(10) unsigned NOT NULL AUTO_INCREMENT;
ALTER TABLE `pham_color` ADD PRIMARY KEY (`ID`);
ALTER TABLE `pham_color` DROP KEY `ID`;
ALTER TABLE `pham_color` CHANGE `PhamID` `Name` int(10) unsigned NOT NULL;
ALTER TABLE `phage` CHANGE `Length` `SequenceLength` mediumint(9) NOT NULL;
ALTER TABLE `phage` CHANGE `Subcluster` `Subcluster2` varchar(5) DEFAULT NULL;
ALTER TABLE `phage` CHANGE `Cluster` `Cluster2` varchar(5) DEFAULT NULL;
ALTER TABLE `phage` ADD COLUMN `Cluster` varchar(5) DEFAULT NULL AFTER `GC`;
UPDATE `phage` SET `Cluster` = `Subcluster2`;
UPDATE `phage` SET `Cluster` = `Cluster2` WHERE `Cluster` is NULL;
ALTER TABLE `phage` CHANGE `HostGenus` `HostStrain` varchar(50) DEFAULT NULL;
UPDATE `version` SET `SchemaVersion` = 6;

### DATA_LOSS_SUMMARY
# LOST_DATA:gene.Parts
