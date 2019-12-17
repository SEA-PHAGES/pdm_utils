# MySQL script to upgrade the database schema from version 3 to 4.
CREATE TABLE `tmrna` (
  `PhageID` varchar(25) NOT NULL,
  `TmrnaID` varchar(35) NOT NULL,
  `LocusTag` varchar(35) DEFAULT NULL,
  `Start` mediumint(9) NOT NULL,
  `Stop` mediumint(9) NOT NULL,
  `Orientation` enum('F','R') NOT NULL,
  `Note` blob,
  `PeptideTag` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`TmrnaID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `trna` (
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
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `trna_structures` (
  `Sequence` varchar(100) NOT NULL,
  `Structure` varchar(300) NOT NULL,
  PRIMARY KEY (`Sequence`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

ALTER TABLE `gene` MODIFY `translation` VARCHAR(5000) DEFAULT NULL;
UPDATE `version` SET `schema_version` = 4;
