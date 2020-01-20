# MySQL statements to change database structure.
# Every conversion step between schema versions is represented by a
# dictionary, with a key for all statements and a key for all data loss info.

# Key names.
statements = "statements"
inaccurate_column = "INACCURATE_COLUMN"
lost_column = "LOST_COLUMN"
lost_table = "LOST_TABLE"
lost_data = "LOST_DATA"
step_summary_dict = "step_summary_dict"

CONVERSION_STEPS = {

    # Upgrade steps
    "upgrade_0_to_1": {
        statements: [
            """CREATE TABLE `version` (
                  `version` int(11) unsigned NOT NULL,
                  PRIMARY KEY (`version`)
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """INSERT INTO `version` (`version`) VALUES ('0');""",
            """ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_1`;""",
            """ALTER TABLE `gene`
                  ADD CONSTRAINT `gene_ibfk_2`
                  FOREIGN KEY (`PhageID`)
                  REFERENCES `phage` (`PhageID`)
                  ON DELETE CASCADE
                  ON UPDATE CASCADE;""",
            """ALTER TABLE `gene_domain` DROP FOREIGN KEY `gene_domain_ibfk_1`;""",
            """ALTER TABLE `gene_domain`
                  ADD CONSTRAINT `gene_domain_ibfk_1`
                  FOREIGN KEY (`GeneID`)
                  REFERENCES `gene` (`GeneID`)
                  ON DELETE CASCADE
                  ON UPDATE CASCADE;""",
            """ALTER TABLE `pham` DROP FOREIGN KEY `pham_ibfk_1`;""",
            """ALTER TABLE `pham`
                  ADD CONSTRAINT `pham_ibfk_1`
                  FOREIGN KEY (`GeneID`)
                  REFERENCES `gene` (`GeneID`)
                  ON DELETE CASCADE
                  ON UPDATE CASCADE;""",
            """ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_1`;""",
            """ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_2`;""",
            """ALTER TABLE `scores_summary`
                  ADD CONSTRAINT `scores_summary_ibfk_1`
                  FOREIGN KEY (`query`)
                  REFERENCES `gene` (`GeneID`)
                  ON DELETE CASCADE
                  ON UPDATE CASCADE,
                  ADD CONSTRAINT `scores_summary_ibfk_2`
                  FOREIGN KEY (`subject`)
                  REFERENCES `gene` (`GeneID`)
                  ON DELETE CASCADE
                  ON UPDATE CASCADE;""",
            """ALTER TABLE `gene` ADD COLUMN `cdd_status` TINYINT(1) NOT NULL DEFAULT '0' AFTER `blast_status`;"""
            ],
        step_summary_dict: {
            inaccurate_column: [
                "version.version",
                "gene.cdd_status"
                ]
            }
        },

    "upgrade_1_to_2": {
        statements: [
            """ALTER TABLE `phage` ADD COLUMN `AnnotationAuthor` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;""",
            """ALTER TABLE `phage` ADD COLUMN `Program` varchar(10) DEFAULT NULL AFTER `status`;""",
            """ALTER TABLE `phage` ADD COLUMN `AnnotationQC` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;""",
            """ALTER TABLE `phage` ADD COLUMN `RetrieveRecord` tinyint(1) NOT NULL DEFAULT '0' AFTER `status`;"""
            ],
        step_summary_dict: {
            inaccurate_column: [
                "phage.AnnotationAuthor",
                "phage.Program",
                "phage.AnnotationQC",
                "phage.RetrieveRecord"
                ]
            }
        },

    "upgrade_2_to_3": {
        statements: [
            """ALTER TABLE `version` ADD COLUMN `schema_version` int(11) unsigned NOT NULL AFTER `version`;""",
            """ALTER TABLE `phage` ADD COLUMN `Subcluster2` varchar(5) DEFAULT NULL AFTER `AnnotationAuthor`;""",
            """ALTER TABLE `phage` ADD COLUMN `Cluster2` varchar(5) DEFAULT NULL AFTER `AnnotationAuthor`;""",
            """ALTER TABLE `phage` ADD COLUMN `TempCluster2` varchar(5) DEFAULT NULL;""",
            """UPDATE `phage` SET `Subcluster2` = `Cluster` WHERE `Cluster` REGEXP '[0-9]$';""",
            """DROP PROCEDURE IF EXISTS split_subcluster;""",
            # Note: below is the syntax if executing directly in mysql:
            # """DELIMITER //
            #     CREATE PROCEDURE split_subcluster()
            #     BEGIN
            #         DECLARE Counter INT DEFAULT 8;
            #         WHILE Counter > -1 DO
            #             UPDATE `phage` SET `TempCluster2` = SUBSTRING(`Subcluster2`, 1, LENGTH(`Subcluster2`)-Counter);
            #             UPDATE `phage` SET `Cluster2` = `TempCluster2` WHERE `TempCluster2` NOT REGEXP '[0-9]$' and `TempCluster2` != '';
            #             SET Counter = Counter - 1;
            #         END WHILE;
            #     END //
            #     DELIMITER;""",
            """CREATE PROCEDURE split_subcluster()
                BEGIN
                    DECLARE Counter INT DEFAULT 8;
                    WHILE Counter > -1 DO
                        UPDATE `phage` SET `TempCluster2` = SUBSTRING(`Subcluster2`, 1, LENGTH(`Subcluster2`)-Counter);
                        UPDATE `phage` SET `Cluster2` = `TempCluster2` WHERE `TempCluster2` NOT REGEXP '[0-9]$' and `TempCluster2` != '';
                        SET Counter = Counter - 1;
                    END WHILE;
                END;""",
            """CALL split_subcluster;""",
            """DROP PROCEDURE split_subcluster;""",
            """UPDATE `phage` SET `Cluster2` = `Cluster` WHERE `Cluster` NOT REGEXP '[0-9]$';""",
            """ALTER TABLE `phage` DROP COLUMN `TempCluster2`;""",
            """ALTER TABLE `phage` DROP COLUMN `Program`;""",
            """ALTER TABLE `gene` ADD COLUMN `LocusTag` varchar(50) DEFAULT NULL AFTER `cdd_status`;""",
            """UPDATE `version` SET `schema_version` = 3;"""
            ],
        step_summary_dict: {
            lost_column: ["phage.Program"],
            inaccurate_column: ["gene.LocusTag"]
            }
        },

    "upgrade_3_to_4": {
        statements: [
            """CREATE TABLE `tmrna` (
                  `PhageID` varchar(25) NOT NULL,
                  `TmrnaID` varchar(35) NOT NULL,
                  `LocusTag` varchar(35) DEFAULT NULL,
                  `Start` mediumint(9) NOT NULL,
                  `Stop` mediumint(9) NOT NULL,
                  `Orientation` enum('F','R') NOT NULL,
                  `Note` blob,
                  `PeptideTag` varchar(50) DEFAULT NULL,
                  PRIMARY KEY (`TmrnaID`)
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """CREATE TABLE `trna` (
                  `PhageID` varchar(25) NOT NULL,
                  `TrnaID` varchar(35) NOT NULL,
                  `LocusTag` varchar(35) DEFAULT NULL,
                  `Start` mediumint(9) NOT NULL,
                  `Stop` mediumint(9) NOT NULL,
                  `Length` mediumint(9) NOT NULL,
                  `Orientation` enum('F','R') NOT NULL,
                  `Sequence` varchar(100) NOT NULL,
                  `Product` blob,
                  `Note` blob,
                  `AminoAcid` enum('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val','Undet','OTHER') NOT NULL,
                  `Anticodon` varchar(4) NOT NULL,
                  `InfernalScore` decimal(4,2) DEFAULT NULL,
                  PRIMARY KEY (`TrnaID`),
                  KEY `PhageID` (`PhageID`),
                  FOREIGN KEY (`PhageID`) REFERENCES `phage` (`PhageID`)
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """CREATE TABLE `trna_structures` (
                  `Sequence` varchar(100) NOT NULL,
                  `Structure` varchar(300) NOT NULL,
                  PRIMARY KEY (`Sequence`)
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """ALTER TABLE `gene` MODIFY `translation` VARCHAR(5000) DEFAULT NULL;""",
            """UPDATE `version` SET `schema_version` = 4;"""
            ],
        step_summary_dict: {}
        },

    "upgrade_4_to_5": {
        statements: [
            """DROP TABLE `node`;""",
            """DROP TABLE `host_range`;""",
            """DROP TABLE `host`;""",
            """DROP TABLE `pham_history`;""",
            """DROP TABLE `pham_old`;""",
            """DROP TABLE `scores_summary`;""",
            """ALTER TABLE `phage` DROP COLUMN `Prophage`;""",
            """ALTER TABLE `phage` DROP COLUMN `Isolated`;""",
            """ALTER TABLE `phage` DROP COLUMN `ProphageOffset`;""",
            """ALTER TABLE `phage` DROP COLUMN `DateLastSearched`;""",
            """ALTER TABLE `phage` DROP COLUMN `AnnotationQC`;""",
            """ALTER TABLE `gene` DROP COLUMN `StartCodon`;""",
            """ALTER TABLE `gene` DROP COLUMN `StopCodon`;""",
            """ALTER TABLE `gene` DROP COLUMN `GC1`;""",
            """ALTER TABLE `gene` DROP COLUMN `GC2`;""",
            """ALTER TABLE `gene` DROP COLUMN `GC3`;""",
            """ALTER TABLE `gene` DROP COLUMN `GC`;""",
            """ALTER TABLE `gene` DROP COLUMN `LeftNeighbor`;""",
            """ALTER TABLE `gene` DROP COLUMN `RightNeighbor`;""",
            """ALTER TABLE `gene` DROP COLUMN `clustalw_status`;""",
            """ALTER TABLE `gene` DROP COLUMN `blast_status`;""",
            """ALTER TABLE `gene` DROP COLUMN `TypeID`;""",
            """ALTER TABLE `pham` DROP COLUMN `orderAdded`;""",
            """ALTER TABLE `phage` MODIFY `status` varchar(7);""",
            """UPDATE `phage` SET `status` = 'unknown' WHERE `status` = 'gbk';""",
            """ALTER TABLE `phage` MODIFY `status` enum('unknown','draft','final');""",
            """UPDATE `version` SET `schema_version` = 5;"""
            ],
        step_summary_dict: {
            lost_table: [
                "node",
                "host_range",
                "host",
                "pham_history",
                "pham_old",
                "scores_summary"
                ],
            lost_column: [
                "phage.Prophage",
                "phage.Isolated",
                "phage.ProphageOffset",
                "phage.DateLastSearched",
                "phage.AnnotationQC",
                "gene.StartCodon",
                "gene.StopCodon",
                "gene.GC1",
                "gene.GC2",
                "gene.GC3",
                "gene.GC",
                "gene.LeftNeighbor",
                "gene.RightNeighbor",
                "gene.clustalw_status",
                "gene.blast_status",
                "gene.TypeID",
                "pham.orderAdded"
                ]
            }
        },

    "upgrade_5_to_6": {
        statements: [
            """ALTER TABLE `domain` CHANGE `id` `ID` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `domain` CHANGE `description` `Description` blob;""",
            """ALTER TABLE `gene_domain` CHANGE `id` `ID` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `gene_domain` CHANGE `expect` `Expect` double unsigned NOT NULL;""",
            """ALTER TABLE `phage` CHANGE `status` `Status` enum('unknown','draft','final') DEFAULT NULL;""",
            """ALTER TABLE `pham` CHANGE `name` `Name` int(10) unsigned DEFAULT NULL;""",
            """ALTER TABLE `pham_color` CHANGE `id` `ID` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `pham_color` CHANGE `name` `Name` int(10) unsigned NOT NULL;""",
            """ALTER TABLE `pham_color` CHANGE `color` `Color` char(7) NOT NULL;""",
            """ALTER TABLE `version` CHANGE `version` `Version` int(11) unsigned NOT NULL;""",
            """ALTER TABLE `gene` CHANGE `translation` `Translation` varchar(5000) DEFAULT NULL;""",
            """ALTER TABLE `gene` CHANGE `id` `ID` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `gene` CHANGE `cdd_status` `DomainStatus` tinyint(1) NOT NULL DEFAULT '0';""",
            """ALTER TABLE `version` CHANGE `schema_version` `SchemaVersion` int(11) unsigned NOT NULL;""",
            """ALTER TABLE `domain` CHANGE `hit_id` `HitID` varchar(25) NOT NULL;""",
            """ALTER TABLE `gene_domain` CHANGE `hit_id` `HitID` varchar(25) NOT NULL;""",
            """ALTER TABLE `gene_domain` CHANGE `query_start` `QueryStart` int(10) unsigned NOT NULL;""",
            """ALTER TABLE `gene_domain` CHANGE `query_end` `QueryEnd` int(10) unsigned NOT NULL;""",
            """UPDATE `version` SET `SchemaVersion` = 6;"""
            ],
        step_summary_dict: {}
        },

    "upgrade_6_to_7": {
        statements: [
            """ALTER TABLE `phage` CHANGE `HostStrain` `HostGenus` varchar(50) DEFAULT NULL;""",
            """ALTER TABLE `phage` DROP COLUMN `Cluster`;""",
            """ALTER TABLE `phage` CHANGE `Cluster2` `Cluster` varchar(5) DEFAULT NULL;""",
            """ALTER TABLE `phage` CHANGE `Subcluster2` `Subcluster` varchar(5) DEFAULT NULL;""",
            """ALTER TABLE `phage` CHANGE `SequenceLength` `Length` mediumint(9) NOT NULL;""",
            """ALTER TABLE `pham_color` CHANGE `Name` `PhamID` int(10) unsigned NOT NULL;""",
            """ALTER TABLE `pham_color` DROP `ID`;""",
            """ALTER TABLE `pham_color` ADD PRIMARY KEY (`PhamID`);""",
            """ALTER TABLE `gene` DROP `ID`;""",
            """ALTER TABLE `gene` ADD COLUMN `Parts` tinyint(1) DEFAULT NULL;""",
            """ALTER TABLE `gene` ADD COLUMN `PhamID` int(10) unsigned DEFAULT NULL;""",
            """UPDATE `gene`
                LEFT JOIN `pham` ON `gene`.`GeneID` = `pham`.`GeneID`
                SET `gene`.`PhamID` = `pham`.`Name`;""",
            """ALTER TABLE `gene`
                ADD FOREIGN KEY (`PhamID`)
                REFERENCES `pham_color` (`PhamID`);""",
            """DROP TABLE `pham`;""",
            """RENAME TABLE `pham_color` TO `pham`;""",
            """UPDATE `version` SET `SchemaVersion` = 7;"""
            ],
        step_summary_dict: {
            inaccurate_column: ["gene.Parts"]
            }
        },

    "upgrade_7_to_8": {
        statements: [
            """ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_3`;""",
            """ALTER TABLE `gene`
                  ADD FOREIGN KEY (`PhamID`)
                  REFERENCES `pham` (`PhamID`)
                  ON DELETE SET NULL;""",
            """UPDATE `version` SET `SchemaVersion` = 8;"""
            ],
        step_summary_dict: {}
        },

    # Downgrade steps
    "downgrade_1_to_0": {
        statements: [
            """ALTER TABLE `gene` DROP COLUMN `cdd_status`;""",
            """ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_1`;""",
            """ALTER TABLE `scores_summary` DROP FOREIGN KEY `scores_summary_ibfk_2`;""",
            """ALTER TABLE `scores_summary`
                  ADD CONSTRAINT `scores_summary_ibfk_1`
                  FOREIGN KEY (`query`)
                  REFERENCES `gene` (`GeneID`),
                  ADD CONSTRAINT `scores_summary_ibfk_2`
                  FOREIGN KEY (`subject`)
                  REFERENCES `gene` (`GeneID`);""",
            """ALTER TABLE `pham` DROP FOREIGN KEY `pham_ibfk_1`;""",
            """ALTER TABLE `pham`
                  ADD CONSTRAINT `pham_ibfk_1`
                  FOREIGN KEY (`GeneID`)
                  REFERENCES `gene` (`GeneID`);""",
            """ALTER TABLE `gene_domain` DROP FOREIGN KEY `gene_domain_ibfk_1`;""",
            """ALTER TABLE `gene_domain`
                  ADD CONSTRAINT `gene_domain_ibfk_1`
                  FOREIGN KEY (`GeneID`)
                  REFERENCES `gene` (`GeneID`);""",
            """ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_2`;""",
            """ALTER TABLE `gene`
                  ADD CONSTRAINT `gene_ibfk_1`
                  FOREIGN KEY (`PhageID`)
                  REFERENCES `phage` (`PhageID`);""",
            """DROP TABLE `version`;"""
            ],
        step_summary_dict: {
            lost_table: ["version"],
            lost_column: ["gene.cdd_status"]
            }
        },

    "downgrade_2_to_1": {
        statements: [
            """ALTER TABLE `phage` DROP COLUMN `RetrieveRecord`;""",
            """ALTER TABLE `phage` DROP COLUMN `AnnotationQC`;""",
            """ALTER TABLE `phage` DROP COLUMN `Program`;""",
            """ALTER TABLE `phage` DROP COLUMN `AnnotationAuthor`;"""
            ],
        step_summary_dict: {
            lost_column: [
                "phage.RetrieveRecord",
                "phage.AnnotationQC",
                "phage.Program",
                "phage.AnnotationAuthor"
                ]
            }
        },

    "downgrade_3_to_2": {
        statements: [
            """ALTER TABLE `phage` ADD COLUMN `Program` varchar(10) DEFAULT NULL AFTER `AnnotationQC`;""",
            """ALTER TABLE `phage` DROP COLUMN `Cluster2`;""",
            """ALTER TABLE `phage` DROP COLUMN `Subcluster2`;""",
            """ALTER TABLE `gene` DROP COLUMN `LocusTag`;""",
            """ALTER TABLE `version` DROP COLUMN `schema_version`;""",
            ],
        step_summary_dict: {
            lost_column: [
                "gene.LocusTag",
                "version.schema_version"
                ],
            inaccurate_column: ["phage.Program"]
            }
        },

    "downgrade_4_to_3": {
        statements: [
            """UPDATE `gene` SET `translation` = NULL WHERE LENGTH(`translation`) > 4000;""",
            """ALTER TABLE `gene` MODIFY `translation` VARCHAR(4000) DEFAULT NULL;""",
            """DROP TABLE `tmrna`;""",
            """DROP TABLE `trna`;""",
            """DROP TABLE `trna_structures`;""",
            """UPDATE `version` SET `schema_version` = 3;"""
            ],
        step_summary_dict: {
            lost_table: [
                "tmrna",
                "trna",
                "trna_structures"
                ],
            lost_data: ["gene.translation"]
            }
        },

    "downgrade_5_to_4": {
        statements: [
            """ALTER TABLE `phage` MODIFY `status` varchar(7) DEFAULT NULL;""",
            """UPDATE `phage` SET `status` = 'gbk' WHERE `status` = 'unknown';""",
            """ALTER TABLE `phage` MODIFY `status` varchar(5) DEFAULT NULL;""",
            """ALTER TABLE `pham` ADD COLUMN `orderAdded` int(5) unsigned NOT NULL AFTER `Name`;""",
            """ALTER TABLE `pham` ADD KEY (`orderAdded`);""",
            """ALTER TABLE `pham` MODIFY `orderAdded` int(5) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `gene` ADD COLUMN `TypeID` varchar(10) DEFAULT NULL AFTER `Name`;""",
            """UPDATE `gene` SET `TypeID` = 'CDS';""",
            """ALTER TABLE `gene` ADD COLUMN `blast_status` enum('avail','pending','stale','done') NOT NULL DEFAULT 'avail' AFTER `id`;""",
            """ALTER TABLE `gene` ADD COLUMN `clustalw_status` enum('avail','pending','stale','done') NOT NULL DEFAULT 'avail' AFTER `id`;""",
            """ALTER TABLE `gene` ADD COLUMN `RightNeighbor` varchar(25) DEFAULT NULL AFTER `Orientation`;""",
            """ALTER TABLE `gene` ADD COLUMN `LeftNeighbor` varchar(25) DEFAULT NULL AFTER `Orientation`;""",
            """ALTER TABLE `gene` ADD COLUMN `GC` float DEFAULT NULL AFTER `Orientation`;""",
            """ALTER TABLE `gene` ADD COLUMN `GC3` float DEFAULT NULL AFTER `Orientation`;""",
            """ALTER TABLE `gene` ADD COLUMN `GC2` float DEFAULT NULL AFTER  `Orientation`;""",
            """ALTER TABLE `gene` ADD COLUMN `GC1` float DEFAULT NULL AFTER `Orientation`;""",
            """ALTER TABLE `gene` ADD COLUMN `StopCodon` enum('TAA','TAG','TGA') DEFAULT NULL AFTER `translation`;""",
            """ALTER TABLE `gene` ADD COLUMN `StartCodon` enum('ATG','GTG','TTG') DEFAULT NULL AFTER `translation`;""",
            """ALTER TABLE `phage` ADD COLUMN `AnnotationQC` tinyint(1) NOT NULL DEFAULT '0' AFTER `RetrieveRecord`;""",
            """ALTER TABLE `phage` ADD COLUMN `DateLastSearched` datetime DEFAULT NULL AFTER `DateLastModified`;""",
            """ALTER TABLE `phage` ADD COLUMN `ProphageOffset` int(11) DEFAULT NULL AFTER `SequenceLength`;""",
            """ALTER TABLE `phage` ADD COLUMN `Prophage` enum('yes','no') DEFAULT NULL AFTER `SequenceLength`;""",
            """ALTER TABLE `phage` ADD COLUMN `Isolated` varchar(100) DEFAULT NULL AFTER `Name`;""",
            """CREATE TABLE `scores_summary` (
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
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """CREATE TABLE `pham_old` (
                  `GeneID` varchar(35) DEFAULT NULL,
                  `name` int(10) unsigned DEFAULT NULL,
                  `orderAdded` int(5) unsigned NOT NULL AUTO_INCREMENT,
                  PRIMARY KEY (`orderAdded`),
                  KEY (`orderAdded`),
                  KEY (`GeneID`),
                  FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON UPDATE CASCADE
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """CREATE TABLE `pham_history` (
                  `name` int(10) unsigned NOT NULL,
                  `parent` int(10) unsigned NOT NULL,
                  `action` enum('join','split') NOT NULL,
                  `datetime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
                  PRIMARY KEY (`name`,`parent`,`action`),
                  KEY (`parent`)
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """CREATE TABLE `host` (
                  `ID` int(11) NOT NULL AUTO_INCREMENT,
                  `Name` varchar(50) NOT NULL,
                  `Accession` varchar(15) DEFAULT NULL,
                  PRIMARY KEY (`ID`)
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """CREATE TABLE `host_range` (
                  `PhageID` varchar(25) NOT NULL,
                  `host_id` int(11) NOT NULL,
                  `plating_efficiency` double unsigned NOT NULL,
                  `below_detection` tinyint(1) NOT NULL,
                  KEY (`PhageID`),
                  KEY (`host_id`),
                  FOREIGN KEY (`PhageID`) REFERENCES phage (`PhageID`),
                  FOREIGN KEY (`host_id`) REFERENCES host (`ID`)
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """CREATE TABLE `node` (
                  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
                  `platform` varchar(15) DEFAULT NULL,
                  `hostname` varchar(50) DEFAULT NULL,
                  PRIMARY KEY (`id`),
                  KEY (`hostname`)
                ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """UPDATE `version` SET `schema_version` = 4;"""
            ],
        step_summary_dict: {
            inaccurate_column: [
                "pham.orderAdded",
                "gene.blast_status",
                "gene.clustalw_status",
                "gene.RightNeighbor",
                "gene.LeftNeighbor",
                "gene.GC",
                "gene.GC3",
                "gene.GC2",
                "gene.GC1",
                "gene.StopCodon",
                "gene.StartCodon",
                "phage.AnnotationQC",
                "phage.DateLastSearched",
                "phage.ProphageOffset",
                "phage.Prophage",
                "phage.Isolated"
                ]
            }
        },

    "downgrade_6_to_5": {
        statements: [
            """ALTER TABLE `gene_domain` CHANGE `QueryEnd` `query_end` int(10) unsigned NOT NULL;""",
            """ALTER TABLE `gene_domain` CHANGE `QueryStart` `query_start` int(10) unsigned NOT NULL;""",
            """ALTER TABLE `gene_domain` CHANGE `HitID` `hit_id` varchar(25) NOT NULL;""",
            """ALTER TABLE `domain` CHANGE `HitID` `hit_id` varchar(25) NOT NULL;""",
            """ALTER TABLE `version` CHANGE `SchemaVersion` `schema_version` int(11) unsigned NOT NULL;""",
            """ALTER TABLE `gene` CHANGE `DomainStatus` `cdd_status` tinyint(1) NOT NULL DEFAULT '0';""",
            """ALTER TABLE `gene` CHANGE `ID` `id` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `gene` CHANGE `Translation` `translation` varchar(5000) DEFAULT NULL;""",
            """ALTER TABLE `version` CHANGE `Version` `version` int(11) unsigned NOT NULL;""",
            """ALTER TABLE `pham_color` CHANGE `Color` `color` char(7) NOT NULL;""",
            """ALTER TABLE `pham_color` CHANGE `Name` `name` int(10) unsigned NOT NULL;""",
            """ALTER TABLE `pham_color` CHANGE `ID` `id` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `pham` CHANGE `Name` `name` int(10) unsigned DEFAULT NULL;""",
            """ALTER TABLE `phage` CHANGE `Status` `status` enum('unknown','draft','final') DEFAULT NULL;""",
            """ALTER TABLE `gene_domain` CHANGE `Expect` `expect` double unsigned NOT NULL;""",
            """ALTER TABLE `gene_domain` CHANGE `ID` `id` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `domain` CHANGE `Description` `description` blob;""",
            """ALTER TABLE `domain` CHANGE `ID` `id` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """UPDATE `version` SET `schema_version` = 5;"""
            ],
        step_summary_dict: {}
        },

    "downgrade_7_to_6": {
        statements: [
            """RENAME TABLE `pham` TO `pham_color`;""",
            """CREATE TABLE `pham` (
                  `GeneID` varchar(35) NOT NULL DEFAULT '',
                  `Name` int(10) unsigned DEFAULT NULL,
                  PRIMARY KEY (`GeneID`),
                  KEY `name_index` (`Name`),
                  CONSTRAINT `pham_ibfk_1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE
                  ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
            """ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_3`;""",
            """INSERT INTO `pham` (`GeneID`, `Name`)
                    SELECT `GeneID`,`PhamID` FROM `gene`;""",
            """ALTER TABLE `gene` DROP `PhamID`;""",
            """ALTER TABLE `gene` DROP `Parts`;""",
            """ALTER TABLE `gene` ADD COLUMN `ID` int(10) unsigned NOT NULL AFTER `Notes`;""",
            """ALTER TABLE `gene` ADD KEY `id` (`ID`);""",
            """ALTER TABLE `gene` MODIFY `ID` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `pham_color` DROP PRIMARY KEY;""",
            """ALTER TABLE `pham_color` ADD COLUMN `ID` int(10) unsigned NOT NULL FIRST;""",
            """ALTER TABLE `pham_color` ADD KEY (`ID`);""",
            """ALTER TABLE `pham_color` MODIFY `ID` int(10) unsigned NOT NULL AUTO_INCREMENT;""",
            """ALTER TABLE `pham_color` ADD PRIMARY KEY (`ID`);""",
            """ALTER TABLE `pham_color` DROP KEY `ID`;""",
            """ALTER TABLE `pham_color` CHANGE `PhamID` `Name` int(10) unsigned NOT NULL;""",
            """ALTER TABLE `phage` CHANGE `Length` `SequenceLength` mediumint(9) NOT NULL;""",
            """ALTER TABLE `phage` CHANGE `Subcluster` `Subcluster2` varchar(5) DEFAULT NULL;""",
            """ALTER TABLE `phage` CHANGE `Cluster` `Cluster2` varchar(5) DEFAULT NULL;""",
            """ALTER TABLE `phage` ADD COLUMN `Cluster` varchar(5) DEFAULT NULL AFTER `GC`;""",
            """UPDATE `phage` SET `Cluster` = `Subcluster2`;""",
            """UPDATE `phage` SET `Cluster` = `Cluster2` WHERE `Cluster` is NULL;""",
            """ALTER TABLE `phage` CHANGE `HostGenus` `HostStrain` varchar(50) DEFAULT NULL;""",
            """UPDATE `version` SET `SchemaVersion` = 6;"""
            ],
        step_summary_dict: {
            lost_data: ["gene.Parts"]
            }
        },

    "downgrade_8_to_7": {
        statements: [
            """ALTER TABLE `gene` DROP FOREIGN KEY `gene_ibfk_3`;""",
            """ALTER TABLE `gene`
                  ADD FOREIGN KEY (`PhamID`)
                  REFERENCES `pham` (`PhamID`);""",
            """UPDATE `version` SET `SchemaVersion` = 7;"""
            ],
        step_summary_dict: {}
        }
    }
