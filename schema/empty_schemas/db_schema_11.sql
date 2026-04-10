-- MySQL dump 10.13  Distrib 8.0.40, for Linux (x86_64)
--
-- Host: localhost    Database: Actino_Draft
-- ------------------------------------------------------
-- Server version	8.0.40-0ubuntu0.22.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `domain`
--

DROP TABLE IF EXISTS `domain`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `domain` (
  `ID` int unsigned NOT NULL AUTO_INCREMENT,
  `HitID` varchar(25) NOT NULL,
  `Description` blob,
  `DomainID` varchar(10) DEFAULT NULL,
  `Name` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE KEY `hit_id` (`HitID`)
) ENGINE=InnoDB AUTO_INCREMENT=1213503 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene`
--

DROP TABLE IF EXISTS `gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `gene` (
  `GeneID` varchar(35) NOT NULL DEFAULT '',
  `PhageID` varchar(25) NOT NULL,
  `Start` mediumint NOT NULL,
  `Stop` mediumint NOT NULL,
  `Length` mediumint NOT NULL,
  `Name` varchar(50) NOT NULL,
  `Translation` blob,
  `Orientation` enum('F','R') DEFAULT NULL,
  `Notes` blob,
  `DomainStatus` tinyint(1) NOT NULL DEFAULT '0',
  `LocusTag` varchar(50) DEFAULT NULL,
  `Parts` tinyint(1) DEFAULT NULL,
  `PhamID` int unsigned DEFAULT NULL,
  `MembraneStatus` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`GeneID`),
  KEY `PhageID` (`PhageID`),
  KEY `PhamID` (`PhamID`),
  CONSTRAINT `gene_ibfk_2` FOREIGN KEY (`PhageID`) REFERENCES `phage` (`PhageID`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `gene_ibfk_3` FOREIGN KEY (`PhamID`) REFERENCES `pham` (`PhamID`) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene_domain`
--

DROP TABLE IF EXISTS `gene_domain`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `gene_domain` (
  `ID` int unsigned NOT NULL AUTO_INCREMENT,
  `GeneID` varchar(35) DEFAULT NULL,
  `HitID` varchar(25) NOT NULL,
  `QueryStart` int unsigned NOT NULL,
  `QueryEnd` int unsigned NOT NULL,
  `Expect` double unsigned NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE KEY `GeneID__hit_id` (`GeneID`,`HitID`),
  KEY `hit_id` (`HitID`),
  CONSTRAINT `gene_domain_ibfk_1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `gene_domain_ibfk_2` FOREIGN KEY (`HitID`) REFERENCES `domain` (`HitID`)
) ENGINE=InnoDB AUTO_INCREMENT=1244241 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene_transmembrane`
--

DROP TABLE IF EXISTS `gene_transmembrane`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `gene_transmembrane` (
  `ID` int unsigned NOT NULL AUTO_INCREMENT,
  `GeneID` varchar(35) DEFAULT NULL,
  `QueryStart` int unsigned NOT NULL,
  `QueryEnd` int unsigned NOT NULL,
  `Type` enum('signal','transmembrane') DEFAULT NULL,
  `Source` enum('deeptmhmm','sosui') DEFAULT NULL,
  PRIMARY KEY (`ID`),
  KEY `GeneID` (`GeneID`),
  CONSTRAINT `fk_gene_transmembrane_gene1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `gene_transmembrane_ibfk_1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`)
) ENGINE=InnoDB AUTO_INCREMENT=200403 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `phage`
--

DROP TABLE IF EXISTS `phage`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `phage` (
  `PhageID` varchar(25) NOT NULL,
  `Accession` varchar(15) NOT NULL,
  `Name` varchar(50) NOT NULL,
  `HostGenus` varchar(50) DEFAULT NULL,
  `Sequence` mediumblob NOT NULL,
  `Length` mediumint NOT NULL,
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

--
-- Table structure for table `pham`
--

DROP TABLE IF EXISTS `pham`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `pham` (
  `PhamID` int unsigned NOT NULL,
  `Color` char(7) NOT NULL,
  PRIMARY KEY (`PhamID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmrna`
--

DROP TABLE IF EXISTS `tmrna`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tmrna` (
  `GeneID` varchar(35) NOT NULL,
  `PhageID` varchar(25) NOT NULL,
  `Start` mediumint NOT NULL,
  `Stop` mediumint NOT NULL,
  `Length` mediumint NOT NULL,
  `Name` varchar(50) NOT NULL,
  `Orientation` enum('F','R') NOT NULL,
  `Note` blob,
  `LocusTag` varchar(35) DEFAULT NULL,
  `PeptideTag` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`GeneID`),
  KEY `PhageID` (`PhageID`),
  CONSTRAINT `tmrna_ibfk_1` FOREIGN KEY (`PhageID`) REFERENCES `phage` (`PhageID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trna`
--

DROP TABLE IF EXISTS `trna`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trna` (
  `GeneID` varchar(35) NOT NULL,
  `PhageID` varchar(25) NOT NULL,
  `Start` mediumint NOT NULL,
  `Stop` mediumint NOT NULL,
  `Length` mediumint NOT NULL,
  `Name` varchar(50) NOT NULL,
  `Orientation` enum('F','R') NOT NULL,
  `Note` blob,
  `LocusTag` varchar(35) DEFAULT NULL,
  `AminoAcid` enum('Ala','Arg','Asn','Asp','Cys','fMet','Gln','Glu','Gly','His','Ile','Ile2','Leu','Lys','Met','Phe','Pro','Pyl','SeC','Ser','Thr','Trp','Tyr','Val','Stop','OTHER') NOT NULL,
  `Anticodon` varchar(4) NOT NULL,
  `Structure` blob,
  `Source` enum('aragorn','trnascan','both') DEFAULT NULL,
  PRIMARY KEY (`GeneID`),
  KEY `PhageID` (`PhageID`),
  CONSTRAINT `trna_ibfk_1` FOREIGN KEY (`PhageID`) REFERENCES `phage` (`PhageID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `version`
--

DROP TABLE IF EXISTS `version`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `version` (
  `Version` int unsigned NOT NULL,
  `SchemaVersion` int unsigned NOT NULL,
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

-- Dump completed on 2025-01-20 20:13:58
