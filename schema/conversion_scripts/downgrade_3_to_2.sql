# MySQL script to downgrade Phamerator database schema from version 3 to 2.
ALTER TABLE `phage` ADD COLUMN `Program` varchar(10) DEFAULT NULL AFTER `AnnotationQC`;
ALTER TABLE `phage` DROP COLUMN `Cluster2`;
ALTER TABLE `phage` DROP COLUMN `Subcluster2`;
ALTER TABLE `gene` DROP COLUMN `LocusTag`;
ALTER TABLE `version` DROP COLUMN `schema_version`;

### DATA_LOSS_SUMMARY
# LOST_COLUMN:gene.LocusTag
# LOST_COLUMN:version.schema_version
# INACCURATE_COLUMN:phage.Program
