# MySQL script to downgrade Phamerator database schema from version 2 to 1.
ALTER TABLE `phage` DROP COLUMN `RetrieveRecord`;
ALTER TABLE `phage` DROP COLUMN `AnnotationQC`;
ALTER TABLE `phage` DROP COLUMN `Program`;
ALTER TABLE `phage` DROP COLUMN `AnnotationAuthor`;

### DATA_LOSS_SUMMARY
# LOST_COLUMN:phage.RetrieveRecord
# LOST_COLUMN:phage.AnnotationQC
# LOST_COLUMN:phage.Program
# LOST_COLUMN:phage.AnnotationAuthor
