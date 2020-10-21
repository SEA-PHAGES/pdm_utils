"""Functions to retrieve phage genome annotation data."""

import re

from sqlalchemy import select

from pdm_utils.functions import basic
from pdm_utils.functions import querying


# ANNOTATION RETRIEVAL
# -----------------------------------------------------------------------------
def get_genes_from_pham(alchemist, pham):
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID
    phamid_obj = gene_obj.c.PhamID

    genes_query = select([geneid_obj]).where(phamid_obj == pham).distinct()

    genes = querying.first_column(alchemist.engine, genes_query)
    return genes


def get_distinct_phams_from_genes(alchemist, geneids):
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID
    phamid_obj = gene_obj.c.PhamID

    phams_query = select([phamid_obj]).distinct()

    phams = querying.first_column(alchemist.engine, phams_query,
                                  in_column=geneid_obj, values=geneids)
    return phams


def get_distinct_annotations_from_genes(alchemist, geneids):
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID
    notes_obj = gene_obj.c.Notes

    notes_query = select([notes_obj]).distinct()

    bytes_annotations = querying.first_column(alchemist.engine, notes_query,
                                              in_column=geneid_obj,
                                              values=geneids)

    annotations = basic.convert_to_decoded(bytes_annotations)
    return annotations


def get_phams_from_genes(alchemist, geneids):
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID
    phamid_obj = gene_obj.c.PhamID

    pham_query = select([phamid_obj])

    phams = []
    for gene in geneids:
        pham = alchemist.engine.execute(pham_query.where(geneid_obj == gene))\
                                                                    .scalar()
        phams.append(pham)

    return phams


def get_annotations_from_genes(alchemist, geneids):
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID
    notes_obj = gene_obj.c.Notes

    notes_query = select([notes_obj])

    annotations = []
    for gene in geneids:
        note = alchemist.engine.execute(notes_query.where(geneid_obj == gene))\
                                                                    .scalar()

        if note is not None:
            note = note.decode("utf-8")
        annotations.append(note)

    return annotations


def get_count_phams_in_genes(alchemist, geneids, incounts=None):
    phams = get_phams_from_genes(alchemist, geneids)

    pham_histogram = {}
    if incounts is not None:
        pham_histogram = incounts

    basic.increment_histogram(phams, pham_histogram)
    return pham_histogram


def get_count_annotations_in_genes(alchemist, geneids, incounts=None):
    annotations = get_annotations_from_genes(alchemist, geneids)

    annotation_counts = {}
    if incounts is not None:
        annotation_counts = incounts

    basic.increment_histogram(annotations, annotation_counts)
    return annotation_counts


def get_count_annotations_in_pham(alchemist, pham, incounts=None):
    genes = get_genes_from_pham(alchemist, pham)
    count_annotations = get_count_annotations_in_genes(alchemist, genes,
                                                       incounts=incounts)
    return count_annotations


def get_relative_gene(alchemist, geneid, pos):
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID

    geneid_format = re.compile("[\w\W]+_CDS_[0-9]+")
    if not re.match(geneid_format, geneid) is None:
        parsed_geneid = re.split("_", geneid)
    else:
        raise ValueError(f"Passed GeneID {geneid} "
                         "is not of the proper GeneID format")

    gene_num = int(parsed_geneid[2])
    rel_gene_pos = gene_num + pos

    rel_geneid = "_".join(parsed_geneid[:2] + [str(rel_gene_pos)])
    geneid_query = select([geneid_obj]).where(geneid_obj == rel_geneid)
    rel_geneid = alchemist.engine.execute(geneid_query).scalar()

    return rel_geneid


def get_adjacent_genes(alchemist, gene):
    left = get_relative_gene(alchemist, gene, -1)
    right = get_relative_gene(alchemist, gene, 1)

    return (left, right)


def get_genes_adjacent_to_pham(alchemist, pham):
    genes = get_genes_from_pham(alchemist, pham)

    left_genes = []
    right_genes = []
    for gene in genes:
        adjacent_genes = get_adjacent_genes(alchemist, gene)

        if not adjacent_genes[0] is None:
            left_genes.append(adjacent_genes[0])
        if not adjacent_genes[1] is None:
            right_genes.append(adjacent_genes[1])

    return (left_genes, right_genes)


def get_distinct_adjacent_phams(alchemist, pham):
    adjacent_genes = get_genes_adjacent_to_pham(alchemist, pham)

    left_phams = get_phams_from_genes(alchemist, adjacent_genes[0])
    right_phams = get_phams_from_genes(alchemist, adjacent_genes[1])

    return (left_phams, right_phams)


def get_count_adjacent_phams_to_pham(alchemist, pham, incounts=None):
    adjacent_genes = get_genes_adjacent_to_pham(alchemist, pham)

    adjacent_phams = ({}, {})
    if incounts is not None:
        adjacent_phams = incounts

    left_in = adjacent_phams[0]
    get_count_phams_in_genes(alchemist, adjacent_genes[0], incounts=left_in)

    right_in = adjacent_phams[1]
    get_count_phams_in_genes(alchemist, adjacent_genes[1], incounts=right_in)

    return adjacent_phams


def get_count_adjacent_annotations_to_pham(alchemist, pham, incounts=None):
    adjacent_genes = get_genes_adjacent_to_pham(alchemist, pham)

    adjacent_annotations = ({}, {})
    if incounts is not None:
        adjacent_annotations = incounts

    left_in = adjacent_annotations[0]
    get_count_annotations_in_genes(alchemist, adjacent_genes[0],
                                   incounts=left_in)

    right_in = adjacent_annotations[1]
    get_count_annotations_in_genes(alchemist, adjacent_genes[1],
                                   incounts=right_in)

    return adjacent_annotations
