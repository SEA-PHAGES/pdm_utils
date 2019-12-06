# MySQL script to upgrade Phamerator database schema from version 4 to 5.
DROP TABLE `node`;
DROP TABLE `host_range`;
DROP TABLE `host`;
DROP TABLE `pham_history`;
DROP TABLE `pham_old`;
DROP TABLE `scores_summary`;
ALTER TABLE `phage` DROP COLUMN `Prophage`;
ALTER TABLE `phage` DROP COLUMN `Isolated`;
ALTER TABLE `phage` DROP COLUMN `ProphageOffset`;
ALTER TABLE `phage` DROP COLUMN `DateLastSearched`;
ALTER TABLE `phage` DROP COLUMN `AnnotationQC`;
ALTER TABLE `gene` DROP COLUMN `StartCodon`;
ALTER TABLE `gene` DROP COLUMN `StopCodon`;
ALTER TABLE `gene` DROP COLUMN `GC1`;
ALTER TABLE `gene` DROP COLUMN `GC2`;
ALTER TABLE `gene` DROP COLUMN `GC3`;
ALTER TABLE `gene` DROP COLUMN `GC`;
ALTER TABLE `gene` DROP COLUMN `LeftNeighbor`;
ALTER TABLE `gene` DROP COLUMN `RightNeighbor`;
ALTER TABLE `gene` DROP COLUMN `clustalw_status`;
ALTER TABLE `gene` DROP COLUMN `blast_status`;
ALTER TABLE `gene` DROP COLUMN `TypeID`;
ALTER TABLE `pham` DROP COLUMN `orderAdded`;
ALTER TABLE `phage` MODIFY `status` varchar(7);
UPDATE `phage` SET `status` = 'unknown' WHERE `status` = 'gbk';
ALTER TABLE `phage` MODIFY `status` enum('unknown','draft','final');
UPDATE `version` SET `schema_version` = 5;

### DATA_LOSS_SUMMARY
# LOST_TABLE:node
# LOST_TABLE:host_range
# LOST_TABLE:host
# LOST_TABLE:pham_history
# LOST_TABLE:pham_old
# LOST_TABLE:scores_summary
# LOST_COLUMN:phage.Prophage
# LOST_COLUMN:phage.Isolated
# LOST_COLUMN:phage.ProphageOffset
# LOST_COLUMN:phage.DateLastSearched
# LOST_COLUMN:phage.AnnotationQC
# LOST_COLUMN:gene.StartCodon
# LOST_COLUMN:gene.StopCodon
# LOST_COLUMN:gene.GC1
# LOST_COLUMN:gene.GC2
# LOST_COLUMN:gene.GC3
# LOST_COLUMN:gene.GC
# LOST_COLUMN:gene.LeftNeighbor
# LOST_COLUMN:gene.RightNeighbor
# LOST_COLUMN:gene.clustalw_status
# LOST_COLUMN:gene.blast_status
# LOST_COLUMN:gene.TypeID
# LOST_COLUMN:pham.orderAdded
