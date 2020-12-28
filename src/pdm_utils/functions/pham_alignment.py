import shlex
from subprocess import (Popen, DEVNULL)

from Bio import SeqIO

from pdm_utils.functions import (fileio, multithread, parallelize)


def run_clustalo(fasta_path, aln_path):
    command = (f"clustalo -i {fasta_path} --infmt=fasta -o {aln_path} "
               "--outfmt=fasta --output-order=tree-order --threads=1 "
               "--seqtype=protein")

    command = shlex.split(command)
    with Popen(command, stdout=DEVNULL, stderr=DEVNULL) as proc:
        proc.wait()

    return aln_path


def get_all_pham_gene_translations(alchemist):
    engine = alchemist.engine

    # Build phage>>cluster lookup table
    cluster_lookup = dict()
    query = "SELECT PhageID, Cluster, Subcluster FROM phage"
    results = engine.execute(query)
    for result in results:
        phageid = result["PhageID"]
        cluster = result["Cluster"]
        subcluster = result["Subcluster"]
        if cluster is None:
            cluster_lookup[phageid] = "Singleton"
        elif subcluster is None:
            cluster_lookup[phageid] = cluster
        else:
            cluster_lookup[phageid] = subcluster

    # Build pham>>translation>>gene lookup table
    phams = dict()
    query = ("SELECT PhamID, LocusTag, Name, Translation, Notes, PhageID "
             "FROM gene")
    results = engine.execute(query)
    for result in results:
        phageid = result["PhageID"]
        locus = result["LocusTag"]
        translation = result["Translation"].decode('utf-8')
        product = result["Notes"].decode('utf-8')
        phamid = result["PhamID"]
        if locus is not None:
            name = locus.split('_')[-1]
        else:
            name = result["Name"]
        geneid = f"{phageid}_{name}"
        if product is not None and product > "":
            geneid += f" ({product})"
        cluster = cluster_lookup[phageid]

        geneid = "".join(["[", cluster, "]", " ", geneid])

        pham_translations = phams.get(phamid, dict())
        gene_ids = pham_translations.get(translation, [])
        gene_ids.append(geneid)
        pham_translations[translation] = gene_ids
        phams[phamid] = pham_translations

    engine.dispose()

    return phams


def dump_pham_out_fastas(working_dir, phams_dict, neglect_singles=True,
                         threads=1, verbose=False):
    pham_fasta_map = dict()

    if verbose:
        print("...Writing pham gene fasta files...")

    work_items = []
    for pham, pham_translations in phams_dict.items():
        filename = f"{pham}_genes.fasta"
        filepath = working_dir.joinpath(filename)

        gs_to_ts = {}
        for translation, gene_ids in pham_translations.items():
            gs_to_ts[gene_ids[0]] = translation

        if (not neglect_singles) or (len(gs_to_ts) > 1):
            pham_fasta_map[pham] = filepath

        work_items.append((gs_to_ts, filepath))

    multithread.multithread(work_items, threads, fileio.write_fasta,
                            verbose=verbose)

    return pham_fasta_map


def align_pham_out_fastas(working_dir, pham_fasta_map, threads=1,
                          verbose=False):
    pham_aln_map = dict()

    if verbose:
        print("...Aligning pham gene fasta files...")

    work_items = []
    for pham, filepath in pham_fasta_map.items():
        aln_name = filepath.with_suffix(".aln").name
        aln_path = working_dir.joinpath(aln_name)

        pham_aln_map[pham] = aln_path

        work_items.append((filepath, aln_path))

    parallelize.parallelize(work_items, threads, run_clustalo, verbose=verbose)

    return pham_aln_map


def reintroduce_fasta_duplicates(filepath, ts_to_gs):
    with filepath.open(mode="r") as filehandle:
        records = []
        for record in SeqIO.parse(filehandle, "fasta"):
            records.append(record)

    gs_to_ts = {}
    for record in records:
        translation = record.seq
        ungapped_translation = translation.ungap(gap="-")

        geneids = ts_to_gs.get(str(ungapped_translation), [])
        for geneid in geneids:
            gs_to_ts[geneid] = str(translation)

    fileio.write_fasta(gs_to_ts, filepath)


def reintroduce_pham_fasta_duplicates(pham_path_map, pham_translations_dict,
                                      threads=1, verbose=False):
    work_items = []
    for pham, filepath in pham_path_map.items():
        ts_to_gs = pham_translations_dict.get(pham, {})

        work_items.append((filepath, ts_to_gs))

    multithread.multithread(work_items, threads, reintroduce_fasta_duplicates,
                            verbose=verbose)
