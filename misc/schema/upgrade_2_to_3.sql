# MySQL script to upgrade Phamerator database schema from version 2 to 3.
ALTER TABLE version ADD COLUMN schema_version int(11) unsigned NOT NULL AFTER version;
ALTER TABLE phage ADD COLUMN Subcluster2 varchar(5) DEFAULT NULL AFTER AnnotationAuthor;
ALTER TABLE phage ADD COLUMN Cluster2 varchar(5) DEFAULT NULL AFTER AnnotationAuthor;
ALTER TABLE phage DROP COLUMN Program;
UPDATE version SET schema_version = 3;
