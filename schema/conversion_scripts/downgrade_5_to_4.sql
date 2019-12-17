# MySQL script to downgrade the database schema from version 5 to 4.
ALTER TABLE `phage` MODIFY `status` varchar(7) DEFAULT NULL;
UPDATE `phage` SET `status` = 'gbk' WHERE `status` = 'unknown';
ALTER TABLE `phage` MODIFY `status` varchar(5) DEFAULT NULL;
ALTER TABLE `pham` ADD COLUMN `orderAdded` int(5) unsigned NOT NULL AFTER `Name`;
ALTER TABLE `pham` ADD KEY (`orderAdded`);
ALTER TABLE `pham` MODIFY `orderAdded` int(5) unsigned NOT NULL AUTO_INCREMENT;
ALTER TABLE `gene` ADD COLUMN `TypeID` varchar(10) DEFAULT NULL AFTER `Name`;
UPDATE `gene` SET `TypeID` = 'CDS';
ALTER TABLE `gene` ADD COLUMN `blast_status` enum('avail','pending','stale','done') NOT NULL DEFAULT 'avail' AFTER `id`;
ALTER TABLE `gene` ADD COLUMN `clustalw_status` enum('avail','pending','stale','done') NOT NULL DEFAULT 'avail' AFTER `id`;
ALTER TABLE `gene` ADD COLUMN `RightNeighbor` varchar(25) DEFAULT NULL AFTER `Orientation`;
ALTER TABLE `gene` ADD COLUMN `LeftNeighbor` varchar(25) DEFAULT NULL AFTER `Orientation`;
ALTER TABLE `gene` ADD COLUMN `GC` float DEFAULT NULL AFTER `Orientation`;
ALTER TABLE `gene` ADD COLUMN `GC3` float DEFAULT NULL AFTER `Orientation`;
ALTER TABLE `gene` ADD COLUMN `GC2` float DEFAULT NULL AFTER  `Orientation`;
ALTER TABLE `gene` ADD COLUMN `GC1` float DEFAULT NULL AFTER `Orientation`;
ALTER TABLE `gene` ADD COLUMN `StopCodon` enum('TAA','TAG','TGA') DEFAULT NULL AFTER `translation`;
ALTER TABLE `gene` ADD COLUMN `StartCodon` enum('ATG','GTG','TTG') DEFAULT NULL AFTER `translation`;
ALTER TABLE `phage` ADD COLUMN `AnnotationQC` tinyint(1) NOT NULL DEFAULT '0' AFTER `RetrieveRecord`;
ALTER TABLE `phage` ADD COLUMN `DateLastSearched` datetime DEFAULT NULL AFTER `DateLastModified`;
ALTER TABLE `phage` ADD COLUMN `ProphageOffset` int(11) DEFAULT NULL AFTER `SequenceLength`;
ALTER TABLE `phage` ADD COLUMN `Prophage` enum('yes','no') DEFAULT NULL AFTER `SequenceLength`;
ALTER TABLE `phage` ADD COLUMN `Isolated` varchar(100) DEFAULT NULL AFTER `Name`;

CREATE TABLE `scores_summary` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `query` varchar(35) DEFAULT NULL,
  `subject` varchar(35) DEFAULT NULL,
  `blast_score` double unsigned DEFAULT NULL,
  `clustalw_score` decimal(5,4) unsigned DEFAULT NULL,
  `blast_bit_score` double unsigned DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY (`query`),
  KEY (`subject`),
  FOREIGN KEY (`query`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE,
  FOREIGN KEY (`subject`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `pham_old` (
  `GeneID` varchar(35) DEFAULT NULL,
  `name` int(10) unsigned DEFAULT NULL,
  `orderAdded` int(5) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`orderAdded`),
  KEY (`orderAdded`),
  KEY (`GeneID`),
  FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `pham_history` (
  `name` int(10) unsigned NOT NULL,
  `parent` int(10) unsigned NOT NULL,
  `action` enum('join','split') NOT NULL,
  `datetime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`name`,`parent`,`action`),
  KEY (`parent`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `host` (
  `ID` int(11) NOT NULL AUTO_INCREMENT,
  `Name` varchar(50) NOT NULL,
  `Accession` varchar(15) DEFAULT NULL,
  PRIMARY KEY (`ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `host_range` (
  `PhageID` varchar(25) NOT NULL,
  `host_id` int(11) NOT NULL,
  `plating_efficiency` double unsigned NOT NULL,
  `below_detection` tinyint(1) NOT NULL,
  KEY (`PhageID`),
  KEY (`host_id`),
  FOREIGN KEY (`PhageID`) REFERENCES phage (`PhageID`),
  FOREIGN KEY (`host_id`) REFERENCES host (`ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `node` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `platform` varchar(15) DEFAULT NULL,
  `hostname` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY (`hostname`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

UPDATE `version` SET `schema_version` = 4;

### DATA_LOSS_SUMMARY
# INACCURATE_COLUMN:pham.orderAdded
# INACCURATE_COLUMN:gene.blast_status
# INACCURATE_COLUMN:gene.clustalw_status
# INACCURATE_COLUMN:gene.RightNeighbor
# INACCURATE_COLUMN:gene.LeftNeighbor
# INACCURATE_COLUMN:gene.GC
# INACCURATE_COLUMN:gene.GC3
# INACCURATE_COLUMN:gene.GC2
# INACCURATE_COLUMN:gene.GC1
# INACCURATE_COLUMN:gene.StopCodon
# INACCURATE_COLUMN:gene.StartCodon
# INACCURATE_COLUMN:phage.AnnotationQC
# INACCURATE_COLUMN:phage.DateLastSearched
# INACCURATE_COLUMN:phage.ProphageOffset
# INACCURATE_COLUMN:phage.Prophage
# INACCURATE_COLUMN:phage.Isolated
