import sqlalchemy
from sqlalchemy.orm import mapper
from sqlalchemy.orm import relationship

from pdm_utils.classes.genome import Genome
from pdm_utils.classes.cds import Cds
from pdm_utils.functions import parsing

#-----------------------------------------------------------------------------
#AUTOMAPPER FUNCTIONS
def get_map(mapper, table):
    """Get SQLAlchemy ORM map object.

    :param mapper: Connected and prepared SQLAlchemy automap base object.
    :type mapper: DeclarativeMeta
    :param table: Case-insensitive table to retrieve a ORM map for.
    :type table: str
    :returns: SQLAlchemy mapped object.
    :rtype: DeclarativeMeta
    """
    table = parsing.translate_table(mapper.metadata, table)
    return mapper.classes[table]

#-----------------------------------------------------------------------------
#CLASSICAL MAPPING FUNCTIONS

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


