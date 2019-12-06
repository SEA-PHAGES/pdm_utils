# MySQL script to downgrade Phamerator database schema from version 1 to 0.
ALTER TABLE `gene` DROP COLUMN `cdd_status`;

ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_1`;
ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_2`;
ALTER TABLE `scores_summary`
  ADD CONSTRAINT `scores_summary_ibfk_1`
  FOREIGN KEY (`query`)
  REFERENCES `gene` (`GeneID`),
  ADD CONSTRAINT `scores_summary_ibfk_2`
  FOREIGN KEY (`subject`)
  REFERENCES `gene` (`GeneID`);

ALTER TABLE `pham` DROP FOREIGN KEY `pham_ibfk_1`;
ALTER TABLE `pham`
  ADD CONSTRAINT `pham_ibfk_1`
  FOREIGN KEY (`GeneID`)
  REFERENCES `gene` (`GeneID`);

ALTER TABLE `gene_domain` DROP FOREIGN KEY `gene_domain_ibfk_1`;
ALTER TABLE `gene_domain`
  ADD CONSTRAINT `gene_domain_ibfk_1`
  FOREIGN KEY (`GeneID`)
  REFERENCES `gene` (`GeneID`);

ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_2`;
ALTER TABLE `gene`
  ADD CONSTRAINT `gene_ibfk_1`
  FOREIGN KEY (`PhageID`)
  REFERENCES `phage` (`PhageID`);

DROP TABLE `version`;

### DATA_LOSS_SUMMARY
# LOST_TABLE:version
# LOST_COLUMN:gene.cdd_status
