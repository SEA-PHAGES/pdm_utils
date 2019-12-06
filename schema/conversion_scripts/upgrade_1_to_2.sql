# MySQL script to upgrade Phamerator database schema from version 1 to 2.
ALTER TABLE `phage` ADD COLUMN `AnnotationAuthor` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;
ALTER TABLE `phage` ADD COLUMN `Program` varchar(10) DEFAULT NULL AFTER `status`;
ALTER TABLE `phage` ADD COLUMN `AnnotationQC` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;
ALTER TABLE `phage` ADD COLUMN `RetrieveRecord` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;

### DATA_LOSS_SUMMARY
# INACCURATE_COLUMN:phage.AnnotationAuthor
# INACCURATE_COLUMN:phage.Program
# INACCURATE_COLUMN:phage.AnnotationQC
# INACCURATE_COLUMN:phage.RetrieveRecord
