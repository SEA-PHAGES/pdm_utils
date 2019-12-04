# MySQL script to upgrade Phamerator database schema from version 1 to 2.
# No loss in data.
# Note: Data in several newly-created columns will not be accurate.
ALTER TABLE `phage` ADD COLUMN `AnnotationAuthor` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;
ALTER TABLE `phage` ADD COLUMN `Program` varchar(10) DEFAULT NULL AFTER `status`;
ALTER TABLE `phage` ADD COLUMN `AnnotationQC` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;
ALTER TABLE `phage` ADD COLUMN `RetrieveRecord` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;
