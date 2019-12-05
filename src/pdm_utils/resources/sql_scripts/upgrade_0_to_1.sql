# MySQL script to upgrade Phamerator database schema from version 0 to 1.
# Note: Data in newly-created gene.cdd_status will not be accurate.
CREATE TABLE `version` (
  `version` int(11) unsigned NOT NULL,
  PRIMARY KEY (`version`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
INSERT INTO `version` (`version`) VALUES ('0');

ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_1`;
ALTER TABLE `gene`
  ADD CONSTRAINT `gene_ibfk_2`
  FOREIGN KEY (`PhageID`)
  REFERENCES `phage` (`PhageID`)
  ON DELETE CASCADE
  ON UPDATE CASCADE;

ALTER TABLE `gene_domain` DROP FOREIGN KEY `gene_domain_ibfk_1`;
ALTER TABLE `gene_domain`
  ADD CONSTRAINT `gene_domain_ibfk_1`
  FOREIGN KEY (`GeneID`)
  REFERENCES `gene` (`GeneID`)
  ON DELETE CASCADE
  ON UPDATE CASCADE;

ALTER TABLE `pham` DROP FOREIGN KEY `pham_ibfk_1`;
ALTER TABLE `pham`
  ADD CONSTRAINT `pham_ibfk_1`
  FOREIGN KEY (`GeneID`)
  REFERENCES `gene` (`GeneID`)
  ON DELETE CASCADE
  ON UPDATE CASCADE;

ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_1`;
ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_2`;
ALTER TABLE `scores_summary`
  ADD CONSTRAINT `scores_summary_ibfk_1`
  FOREIGN KEY (`query`)
  REFERENCES `gene` (`GeneID`)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  ADD CONSTRAINT `scores_summary_ibfk_2`
  FOREIGN KEY (`subject`)
  REFERENCES `gene` (`GeneID`)
  ON DELETE CASCADE
  ON UPDATE CASCADE;

ALTER TABLE `gene` ADD COLUMN `cdd_status` TINYINT(1) NOT NULL DEFAULT '0' AFTER `blast_status`;
