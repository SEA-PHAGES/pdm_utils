# MySQL script to downgrade Phamerator database schema from version 4 to 3.
# Note: Data in several tables will be lost.
DROP TABLE tmrna;
DROP TABLE trna;
DROP TABLE trna_structures;
UPDATE version SET schema_version = 3;
