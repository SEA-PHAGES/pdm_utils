# MySQL script to upgrade Phamerator database schema from version 4 to 5.
# Note: Data in several columns and tables will be lost.
DROP TABLE node;
DROP TABLE host_range;
DROP TABLE host;
DROP TABLE pham_history;
DROP TABLE pham_old;
DROP TABLE scores_summary;
ALTER TABLE phage DROP COLUMN Prophage;
ALTER TABLE phage DROP COLUMN Isolated;
ALTER TABLE phage DROP COLUMN ProphageOffset;
ALTER TABLE phage DROP COLUMN DateLastSearched;
ALTER TABLE phage DROP COLUMN AnnotationQC;
ALTER TABLE gene DROP COLUMN StartCodon;
ALTER TABLE gene DROP COLUMN StopCodon;
ALTER TABLE gene DROP COLUMN GC1;
ALTER TABLE gene DROP COLUMN GC2;
ALTER TABLE gene DROP COLUMN GC3;
ALTER TABLE gene DROP COLUMN GC;
ALTER TABLE gene DROP COLUMN LeftNeighbor;
ALTER TABLE gene DROP COLUMN RightNeighbor;
ALTER TABLE gene DROP COLUMN clustalw_status;
ALTER TABLE gene DROP COLUMN blast_status;
ALTER TABLE gene DROP COLUMN TypeID;
ALTER TABLE pham DROP COLUMN orderAdded;
ALTER TABLE phage MODIFY status varchar(7);
UPDATE phage SET status = 'unknown' WHERE status = 'gbk';
ALTER TABLE phage MODIFY status enum('unknown','draft','final');
UPDATE version SET schema_version = 5;
