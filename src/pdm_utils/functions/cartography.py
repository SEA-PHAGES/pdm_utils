import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String, Text, Table
from sqlalchemy.types import LargeBinary
from sqlalchemy import MetaData, select, ForeignKey, join
from sqlalchemy.orm import sessionmaker, mapper, column_property, relationship
from sqlalchemy.orm import reconstructor
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql import func
from pdm_utils.classes.genome import Genome
from pdm_utils.classes.cds import Cds
from Bio.Seq import Seq
    
def map_cds(metadata):
    gene_table = metadata.tables["gene"]

    cds_properties = {"id" : gene_table.c.GeneID,
                      "genome_id" : gene_table.c.PhageID,
                      "start" : gene_table.c.Start,
                      "stop"  : gene_table.c.Stop,
                      "name"  : gene_table.c.Name,
                      "translation": gene_table.c.Translation,
                      "orientation" : gene_table.c.Orientation,
                      "description" : gene_table.c.Notes,
                      "locus_tag"   : gene_table.c.LocusTag}

    cds_exclusions = [gene_table.c.DomainStatus, gene_table.c.ID]

    map = mapper(Cds, gene_table, properties=cds_properties,
                 exclude_properties=cds_exclusions)

    return map

def map_genome(metadata):
    phage_table = metadata.tables["phage"]
    map_cds(metadata)
   
    genome_properties = {"id"                 : phage_table.c.PhageID,
                         "accession"          : phage_table.c.Accession,
                         "name"               : phage_table.c.Name,
                         "host_genus"         : phage_table.c.HostStrain,
                         "sequence"           : phage_table.c.Sequence,
                         "length"             : phage_table.c.SequenceLength,
                         "date"               : phage_table.c.DateLastModified,
                         "gc"                 : phage_table.c.GC,
                         "annotation_status"  : phage_table.c.Status,
                         "retrieve_record"    : phage_table.c.RetrieveRecord,
                         "annotation_author"  : phage_table.c.AnnotationAuthor,
                         "cluster"            : phage_table.c.Cluster2,
                         "subcluster"         : phage_table.c.Subcluster2,
                         "cds_features"       : relationship(Cds)}

    genome_exclusions = [phage_table.c.Cluster, phage_table.c.Notes]

    map = mapper(Genome, phage_table, properties=genome_properties,
                 exclude_properties=genome_exclusions)

    return map

MAPPINGS = {"genome" : map_genome,
            "cds"    : map_cds}

MAP_TEMPLATES = MAPPINGS.keys()

def get_map(metadata, template):
    mapper = MAPPINGS[template]
    
    map = mapper(metadata)
    return map

