-- MySQL dump 10.13  Distrib 5.5.43, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: travis
-- ------------------------------------------------------
-- Server version	5.5.43-0ubuntu0.14.04.1

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

--
-- Table structure for table `domain`
--

DROP TABLE IF EXISTS `domain`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `domain` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `hit_id` varchar(25) NOT NULL,
  `description` blob,
  `DomainID` varchar(10) DEFAULT NULL,
  `Name` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `hit_id` (`hit_id`)
) ENGINE=InnoDB AUTO_INCREMENT=813827 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene`
--

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
  `TypeID` varchar(10) DEFAULT NULL,
  `translation` varchar(4000) DEFAULT NULL,
  `StartCodon` enum('ATG','GTG','TTG') DEFAULT NULL,
  `StopCodon` enum('TAA','TAG','TGA') DEFAULT NULL,
  `Orientation` enum('F','R') DEFAULT NULL,
  `GC1` float DEFAULT NULL,
  `GC2` float DEFAULT NULL,
  `GC3` float DEFAULT NULL,
  `GC` float DEFAULT NULL,
  `LeftNeighbor` varchar(25) DEFAULT NULL,
  `RightNeighbor` varchar(25) DEFAULT NULL,
  `Notes` blob,
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `clustalw_status` enum('avail','pending','stale','done') NOT NULL DEFAULT 'avail',
  `blast_status` enum('avail','pending','stale','done') NOT NULL DEFAULT 'avail',
  `cdd_status` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`GeneID`),
  KEY `PhageID` (`PhageID`),
  KEY `id` (`id`),
  CONSTRAINT `gene_ibfk_2` FOREIGN KEY (`PhageID`) REFERENCES `phage` (`PhageID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=135628 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene_domain`
--

DROP TABLE IF EXISTS `gene_domain`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene_domain` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `GeneID` varchar(35) DEFAULT NULL,
  `hit_id` varchar(25) NOT NULL,
  `query_start` int(10) unsigned NOT NULL,
  `query_end` int(10) unsigned NOT NULL,
  `expect` double unsigned NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `GeneID__hit_id` (`GeneID`,`hit_id`),
  KEY `hit_id` (`hit_id`),
  CONSTRAINT `gene_domain_ibfk_1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `gene_domain_ibfk_2` FOREIGN KEY (`hit_id`) REFERENCES `domain` (`hit_id`)
) ENGINE=InnoDB AUTO_INCREMENT=980216 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `host`
--

DROP TABLE IF EXISTS `host`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `host` (
  `ID` int(11) NOT NULL AUTO_INCREMENT,
  `Name` varchar(50) NOT NULL,
  `Accession` varchar(15) DEFAULT NULL,
  PRIMARY KEY (`ID`)
) ENGINE=InnoDB AUTO_INCREMENT=4 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `host_range`
--

DROP TABLE IF EXISTS `host_range`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `host_range` (
  `PhageID` varchar(25) NOT NULL,
  `host_id` int(11) NOT NULL,
  `plating_efficiency` double unsigned NOT NULL,
  `below_detection` tinyint(1) NOT NULL,
  KEY `PhageID` (`PhageID`),
  KEY `host_id` (`host_id`),
  CONSTRAINT `host_range_ibfk_1` FOREIGN KEY (`PhageID`) REFERENCES `phage` (`PhageID`),
  CONSTRAINT `host_range_ibfk_2` FOREIGN KEY (`host_id`) REFERENCES `host` (`ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `node`
--

DROP TABLE IF EXISTS `node`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `node` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `platform` varchar(15) DEFAULT NULL,
  `hostname` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `hostname_index` (`hostname`)
) ENGINE=InnoDB AUTO_INCREMENT=2536 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `phage`
--

DROP TABLE IF EXISTS `phage`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phage` (
  `PhageID` varchar(25) NOT NULL,
  `Accession` varchar(15) NOT NULL,
  `Name` varchar(50) NOT NULL,
  `Isolated` varchar(100) DEFAULT NULL,
  `HostStrain` varchar(50) DEFAULT NULL,
  `Sequence` mediumblob NOT NULL,
  `SequenceLength` mediumint(9) NOT NULL,
  `Prophage` enum('yes','no') DEFAULT NULL,
  `ProphageOffset` int(11) DEFAULT NULL,
  `DateLastModified` datetime DEFAULT NULL,
  `DateLastSearched` datetime DEFAULT NULL,
  `Notes` blob,
  `GC` float DEFAULT NULL,
  `Cluster` varchar(5) DEFAULT NULL,
  `status` varchar(5) DEFAULT NULL,
  PRIMARY KEY (`PhageID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pham`
--

DROP TABLE IF EXISTS `pham`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pham` (
  `GeneID` varchar(35) NOT NULL DEFAULT '',
  `name` int(10) unsigned DEFAULT NULL,
  `orderAdded` int(5) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`GeneID`),
  KEY `orderAdded_index` (`orderAdded`),
  KEY `name_index` (`name`),
  CONSTRAINT `pham_ibfk_1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=93667 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pham_color`
--

DROP TABLE IF EXISTS `pham_color`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pham_color` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` int(10) unsigned NOT NULL,
  `color` char(7) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=5677 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pham_history`
--

DROP TABLE IF EXISTS `pham_history`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pham_history` (
  `name` int(10) unsigned NOT NULL,
  `parent` int(10) unsigned NOT NULL,
  `action` enum('join','split') NOT NULL,
  `datetime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`name`,`parent`,`action`),
  KEY `parent_index` (`parent`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pham_old`
--

DROP TABLE IF EXISTS `pham_old`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pham_old` (
  `GeneID` varchar(35) DEFAULT NULL,
  `name` int(10) unsigned DEFAULT NULL,
  `orderAdded` int(5) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`orderAdded`),
  KEY `orderAdded_index` (`orderAdded`),
  KEY `GeneID_index` (`GeneID`),
  CONSTRAINT `pham_old_ibfk_1` FOREIGN KEY (`GeneID`) REFERENCES `gene` (`GeneID`) ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `scores_summary`
--

DROP TABLE IF EXISTS `scores_summary`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `scores_summary` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `query` varchar(35) DEFAULT NULL,
  `subject` varchar(35) DEFAULT NULL,
  `blast_score` double unsigned DEFAULT NULL,
  `clustalw_score` decimal(5,4) unsigned DEFAULT NULL,
  `blast_bit_score` double unsigned DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `scores_summary_ibfk_1` (`query`),
  KEY `scores_summary_ibfk_2` (`subject`)
) ENGINE=InnoDB AUTO_INCREMENT=6106 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `version`
--

DROP TABLE IF EXISTS `version`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `version` (
  `version` int(11) unsigned NOT NULL,
  PRIMARY KEY (`version`)
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

-- Dump completed on 2015-05-07 14:41:38
