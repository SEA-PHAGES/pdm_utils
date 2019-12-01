# MySQL script to downgrade Phamerator database schema from version 7 to 6.
RENAME TABLE pham TO pham_color;
RENAME TABLE cds_pham TO pham;
RENAME TABLE cds_domain TO gene_domain;
RENAME TABLE cds TO gene;
ALTER TABLE phage CHANGE Subcluster Subcluster2  varchar(5) DEFAULT NULL;
ALTER TABLE phage CHANGE Cluster Cluster2 varchar(5) DEFAULT NULL;
ALTER TABLE phage ADD COLUMN Cluster varchar(5) DEFAULT NULL AFTER GC;
UPDATE phage SET Cluster = Subcluster2;
UPDATE phage SET Cluster = Cluster2 WHERE Cluster is NULL;
ALTER TABLE phage CHANGE HostGenus HostStrain varchar(50) DEFAULT NULL;
UPDATE version SET SchemaVersion = 6;
