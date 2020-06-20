"""Functions to retrieve phage genome annotation data."""

import re

import sqlalchemy
from sqlalchemy import and_
from sqlalchemy import select
from sqlalchemy.sql import func

from pdm_utils.functions import querying

def get_relative_gene(alchemist, geneid, pos):
    gene_obj = alchemist.metadata.tables["gene"] 

    geneid_obj = gene_obj.c.GeneID

    geneid_format = re.compile("\w+_CDS_[0-9]+")
    if not re.match(geneid_format, geneid) is None:
        parsed_geneid = re.split("_", geneid)
    else:
        raise ValueError("Passed GeneID is not of the proper GeneID format")

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

def get_genes_from_pham(alchemist, pham):
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID
    phamid_obj = gene_obj.c.PhamID

    genes_query = select([geneid_obj]).where(phamid_obj == pham).distinct()
    
    genes = querying.first_column(alchemist.engine, genes_query)
    return genes

def get_phams_from_genes(alchemist, geneids):
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID
    phamid_obj = gene_obj.c.PhamID

    phams_query = select([phamid_obj]).distinct()
 
    phams = querying.first_column(alchemist.engine, phams_query,
                                  in_column=geneid_obj, values=geneids)

    return phams

def get_genes_adjacent_to_pham(alchemist, pham):
    genes = get_genes_from_pham(alchemist, pham)
    
    left_genes = []
    right_genes = []
    for gene in genes:
        try:
            adjacent_genes = get_adjacent_genes(alchemist, gene)
        except:
            continue

        left = adjacent_genes[0]
        right = adjacent_genes[1]

        if not left is None:
            left_genes.append(left)
        if not right is None:
            right_genes.append(right)

    return (left_genes, right_genes)

def get_adjacent_phams(alchemist, pham):
    adjacent_genes = get_genes_adjacent_to_pham(alchemist, pham)     

    left_phams = get_phams_from_genes(alchemist, adjacent_genes[0])
    right_phams = get_phams_from_genes(alchemist, adjacent_genes[1])

    return (left_phams, right_phams)

def get_count_adjacent_phams(alchemist, pham, incounts=None):
    adjacent_genes = get_genes_adjacent_to_pham(alchemist, pham) 

    count_adjacent_phams_dicts = {}
    if not incounts is None:
        count_adjacent_phams_dicts = incounts

    left_phams_in = None
    if not incounts is None:
        left_phams_in = incounts['left']
    count_adjacent_phams_dicts['left'] = get_count_phams_in_genes(
                            alchemist, left_genes, incounts=left_phams_in)
   
    right_phams_in = None
    if not incounts is None:
        right_phams_in = incounts['right']
    count_adjacent_phams_dicts['right'] = get_count_phams_in_genes(
                            alchemist, right_genes, incounts=right_phams_in)

    return count_adjacent_phams_dicts
    
def get_count_phams_in_genes(alchemist, geneids, incounts=None): 
    gene_obj = alchemist.metadata.tables["gene"]

    geneid_obj = gene_obj.c.GeneID
    phamid_obj = gene_obj.c.PhamID

    pham_query = select([phamid_obj])

    pham_histogram = {}
    if not incounts is None:
        pham_histogram = incounts

    for gene in geneids:
        pham = alchemist.engine.execute(pham_query.where(geneid_obj == gene))\
                                                                    .scalar()
        try:
            pham_histogram[pham] += 1
        except:
            pham_histogram[pham] = 1

    return pham_histogram

def get_count_pham_annotations(alchemist, pham, incounts=None):
    gene_obj = alchemist.metadata.tables["gene"]  

    geneid_obj = gene_obj.c.GeneID
    phamid_obj = gene_obj.c.PhamID
    notes_obj = gene_obj.c.Notes

    annotation_query = select([notes_obj]).where(phamid_obj == pham).distinct()
    annotations = querying.first_column(alchemist.engine, annotation_query)

    annotation_counts = {}
    if not incounts is None:
        annotation_counts = incounts

    for byte_annotation in annotations:
        annotation = None
        if not byte_annotation is None:
            annotation = byte_annotation.decode("utf-8")

        count_query = select([func.count(geneid_obj)]).where(and_(*[
                        (notes_obj == byte_annotation),(phamid_obj == pham)]))
        count = alchemist.engine.execute(count_query).scalar()

        try:
            annotation_counts[annotation] += count
        except:
            annotation_counts[annotation] = count

    return annotation_counts
