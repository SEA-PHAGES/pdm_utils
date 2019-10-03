#This statement adds the version table to the database for version tracking.
CREATE TABLE `version` (
 `version` int(11) unsigned NOT NULL,
 PRIMARY KEY (`version`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

#Execute these statements in SQL to properly update your sql schema for easy genome removal.

ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_1` ;
ALTER TABLE `gene` 
  ADD CONSTRAINT `gene_ibfk_1`
  FOREIGN KEY (`PhageID` )
  REFERENCES `phage` (`PhageID` )
  ON DELETE CASCADE
  ON UPDATE CASCADE;

ALTER TABLE `gene_domain` DROP FOREIGN KEY `gene_domain_ibfk_1` ;
ALTER TABLE `gene_domain` 
  ADD CONSTRAINT `gene_domain_ibfk_1`
  FOREIGN KEY (`GeneID` )
  REFERENCES `gene` (`GeneID` )
  ON DELETE CASCADE
  ON UPDATE CASCADE;

#NOTE!  Some database for some reason have this key as pham_ibfk_2 - change accordingly
ALTER TABLE `pham` DROP FOREIGN KEY `pham_ibfk_1` ;
ALTER TABLE `pham` 
  ADD CONSTRAINT `pham_ibfk_1`
  FOREIGN KEY (`GeneID` )
  REFERENCES `gene` (`GeneID` )
  ON DELETE CASCADE
  ON UPDATE CASCADE;

ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_1` , DROP FOREIGN KEY `scores_summary_ibfk_2` ;
ALTER TABLE `scores_summary` 
  ADD CONSTRAINT `scores_summary_ibfk_1`
  FOREIGN KEY (`query` )
  REFERENCES `gene` (`GeneID` )
  ON DELETE CASCADE
  ON UPDATE CASCADE, 
  ADD CONSTRAINT `scores_summary_ibfk_2`
  FOREIGN KEY (`subject` )
  REFERENCES `gene` (`GeneID` )
  ON DELETE CASCADE
  ON UPDATE CASCADE;

#This statement shrinks database size substantially, as these alignments are no longer used for phameration

TRUNCATE TABLE `scores_summary`;

#This statement adds the framework necessary to use cdd_pp.py

ALTER TABLE `gene` ADD COLUMN `cdd_status` TINYINT(1) NOT NULL  AFTER `blast_status` ;

#This statement will populate the cdd_status column depending on what you need.

SET SQL_SAFE_UPDATES = 0;
update gene
set cdd_status = 0;
 #Or 1 if already completed for current DB
SET SQL_SAFE_UPDATES = 1;

# 01/11/2019 Schema 4 - we now have multiple genomes each with a gene encoding a protein longer than 4000 residues. Graham approved 
# increasing the max protein size to 5000 until the next time this needs to increase

ALTER TABLE gene MODIFY translation VARCHAR(5000);
COMMIT;

# END FILE

