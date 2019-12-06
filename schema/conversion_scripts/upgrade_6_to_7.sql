# MySQL script to upgrade Phamerator database schema from version 6 to 7.
# No loss in data.
ALTER TABLE `phage` CHANGE `HostStrain` `HostGenus` varchar(50) DEFAULT NULL;
ALTER TABLE `phage` DROP COLUMN `Cluster`;
ALTER TABLE `phage` CHANGE `Cluster2` `Cluster` varchar(5) DEFAULT NULL;
ALTER TABLE `phage` CHANGE `Subcluster2` `Subcluster` varchar(5) DEFAULT NULL;
RENAME TABLE `gene` TO `cds`;
RENAME TABLE `gene_domain` TO `cds_domain`;
RENAME TABLE `pham` TO `cds_pham`;
RENAME TABLE `pham_color` TO `pham`;
UPDATE `version` SET `SchemaVersion` = 7;
