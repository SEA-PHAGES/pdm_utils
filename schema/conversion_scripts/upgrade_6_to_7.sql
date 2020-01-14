# MySQL script to upgrade the database schema from version 6 to 7.
ALTER TABLE `phage` CHANGE `HostStrain` `HostGenus` varchar(50) DEFAULT NULL;
ALTER TABLE `phage` DROP COLUMN `Cluster`;
ALTER TABLE `phage` CHANGE `Cluster2` `Cluster` varchar(5) DEFAULT NULL;
ALTER TABLE `phage` CHANGE `Subcluster2` `Subcluster` varchar(5) DEFAULT NULL;
ALTER TABLE `phage` CHANGE `SequenceLength` `Length` mediumint(9) NOT NULL;
ALTER TABLE `pham_color` CHANGE `Name` `PhamID` int(10) unsigned NOT NULL;
ALTER TABLE `pham_color` DROP `ID`;
ALTER TABLE `pham_color` ADD PRIMARY KEY (`PhamID`);
ALTER TABLE `gene` DROP `ID`;
ALTER TABLE `gene` ADD COLUMN `Parts` tinyint(1) DEFAULT NULL;
ALTER TABLE `gene` ADD COLUMN `PhamID` int(10) unsigned DEFAULT NULL;
UPDATE `gene`
  LEFT JOIN `pham` ON `gene`.`GeneID` = `pham`.`GeneID`
  SET `gene`.`PhamID` = `pham`.`Name`;
ALTER TABLE `gene`
  ADD FOREIGN KEY (`PhamID`)
  REFERENCES `pham_color` (`PhamID`);
DROP TABLE `pham`;
RENAME TABLE `pham_color` TO `pham`;
UPDATE `version` SET `SchemaVersion` = 7;

### DATA_LOSS_SUMMARY
# INACCURATE_COLUMN:gene.Parts
