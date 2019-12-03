# MySQL script to downgrade Phamerator database schema from version 6 to 5.
# No loss in data.
ALTER TABLE gene_domain CHANGE QueryEnd query_end int(10) unsigned NOT NULL;
ALTER TABLE gene_domain CHANGE QueryStart query_start int(10) unsigned NOT NULL;
ALTER TABLE gene_domain CHANGE HitID hit_id varchar(25) NOT NULL;
ALTER TABLE domain CHANGE HitID hit_id varchar(25) NOT NULL;
ALTER TABLE version CHANGE SchemaVersion schema_version int(11) unsigned NOT NULL;
ALTER TABLE gene CHANGE DomainStatus cdd_status tinyint(1) NOT NULL DEFAULT '0';
ALTER TABLE gene CHANGE ID id int(10) unsigned NOT NULL AUTO_INCREMENT;
ALTER TABLE gene CHANGE Translation translation varchar(5000) DEFAULT NULL;
ALTER TABLE version CHANGE Version version int(11) unsigned NOT NULL;
ALTER TABLE pham_color CHANGE Color color char(7) NOT NULL;
ALTER TABLE pham_color CHANGE Name name int(10) unsigned NOT NULL;
ALTER TABLE pham_color CHANGE ID id int(10) unsigned NOT NULL AUTO_INCREMENT;
ALTER TABLE pham CHANGE Name name int(10) unsigned DEFAULT NULL;
ALTER TABLE phage CHANGE Status status enum('unknown','draft','final') DEFAULT NULL;
ALTER TABLE gene_domain CHANGE Expect expect double unsigned NOT NULL;
ALTER TABLE gene_domain CHANGE ID id int(10) unsigned NOT NULL AUTO_INCREMENT;
ALTER TABLE domain CHANGE Description description blob;
ALTER TABLE domain CHANGE ID id int(10) unsigned NOT NULL AUTO_INCREMENT;
UPDATE version SET schema_version = 5;
