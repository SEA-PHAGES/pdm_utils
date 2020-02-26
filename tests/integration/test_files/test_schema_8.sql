
/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;
DROP TABLE IF EXISTS `domain`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `domain` (
  `ID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `HitID` varchar(25) NOT NULL,
  `Description` blob,
  `DomainID` varchar(10) DEFAULT NULL,
  `Name` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE KEY `hit_id` (`HitID`)
) ENGINE=InnoDB AUTO_INCREMENT=2202795 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene` (
  `GeneID` varchar(35) NOT NULL DEFAULT '',
  `PhageID` varchar(25) NOT NULL,
  `Start` mediumint(9) NOT NULL,
  `Stop` mediumint(9) NOT NULL,
  `Length` mediumint(9) NOT NULL,
  `Name` varchar(50) NOT NULL,
  `Translation` varchar(5000) DEFAULT NULL,
  `Orientation` enum('F','R') DEFAULT NULL,
  `Notes` blob,
  `DomainStatus` tinyint(1) NOT NULL DEFAULT '0',
  `LocusTag` varchar(50) DEFAULT NULL,
  `Parts` tinyint(1) DEFAULT NULL,
  `PhamID` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`GeneID`),
  KEY `PhageID` (`PhageID`),
  KEY `PhamID` (`PhamID`),
  CONSTRAINT `gene_ibfk_2` FOREIGN KEY (`PhageID`) REFERENCES `phage` (`PhageID`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `gene_ibfk_3` FOREIGN KEY (`PhamID`) REFERENCES `pham` (`PhamID`) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `gene_domain`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene_domain` (
  `ID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `GeneID` varchar(35) DEFAULT NULL,
  `HitID` varchar(25) NOT NULL,
  `QueryStart` int(10) unsigned NOT NULL,
  `QueryEnd` int(10) unsigned NOT NULL,
  `Expect` double unsigned NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE KEY `GeneID__hit_id` (`GeneID`,`HitID`),
  KEY `hit_id` (`HitID`),
  CONSTRAINT `gene_domain_ibfk_1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `gene_domain_ibfk_2` FOREIGN KEY (`HitID`) REFERENCES `domain` (`HitID`)
) ENGINE=InnoDB AUTO_INCREMENT=1413975 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `phage`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phage` (
  `PhageID` varchar(25) NOT NULL,
  `Accession` varchar(15) NOT NULL,
  `Name` varchar(50) NOT NULL,
  `HostGenus` varchar(50) DEFAULT NULL,
  `Sequence` mediumblob NOT NULL,
  `Length` mediumint(9) NOT NULL,
  `DateLastModified` datetime DEFAULT NULL,
  `Notes` blob,
  `GC` float DEFAULT NULL,
  `Status` enum('unknown','draft','final') DEFAULT NULL,
  `RetrieveRecord` tinyint(1) NOT NULL DEFAULT '0',
  `AnnotationAuthor` tinyint(1) NOT NULL DEFAULT '0',
  `Cluster` varchar(5) DEFAULT NULL,
  `Subcluster` varchar(5) DEFAULT NULL,
  PRIMARY KEY (`PhageID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `pham`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pham` (
  `PhamID` int(10) unsigned NOT NULL,
  `Color` char(7) NOT NULL,
  PRIMARY KEY (`PhamID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `tmrna`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
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
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `trna`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
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
  CONSTRAINT `trna_ibfk_1` FOREIGN KEY (`PhageID`) REFERENCES `phage` (`PhageID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `trna_structures`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trna_structures` (
  `Sequence` varchar(100) NOT NULL,
  `Structure` varchar(300) NOT NULL,
  PRIMARY KEY (`Sequence`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `version`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `version` (
  `Version` int(11) unsigned NOT NULL,
  `SchemaVersion` int(11) unsigned NOT NULL,
  PRIMARY KEY (`Version`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

