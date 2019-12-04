# MySQL script to downgrade Phamerator database schema from version 3 to 2.
# Note: Data in several columns will be lost.
# Note: Data in a newly-created column will not be accurate.
ALTER TABLE `phage` ADD COLUMN `Program` varchar(10) DEFAULT NULL AFTER `AnnotationQC`;
ALTER TABLE `phage` DROP COLUMN `Cluster2`;
ALTER TABLE `phage` DROP COLUMN `Subcluster2`;
ALTER TABLE `version` DROP COLUMN `schema_version`;
