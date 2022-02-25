import shlex
from subprocess import (Popen, DEVNULL)

from pdm_utils.functions import (fileio, multithread, multiprocess)


THREAD_LOCK_TIMEOUT = 10
THREAD_JOIN_TIMEOUT = 30


def run_clustalo(fasta_path, aln_path):
    """Runs Clustal Omega to generate a fasta-formatted  multiple sequence
    alignment file

    :param fasta_path: Path to a fasta file containing sequences to be aligned
    :type fasta_path: Path
    :type fasta_path: str
    :param aln_path: The desired path to the aligned sequences file
    :type aln_path: Path
    :type aln_path: str
    :return: Returns the path to the aligned sequences file
    :rtype: Path
    :rtype: str
    """
    command = (f"clustalo -i {fasta_path} --infmt=fasta -o {aln_path} "
               "--outfmt=fasta --output-order=tree-order --threads=1 "
               "--seqtype=protein")

    command = shlex.split(command)
    with Popen(command, stdout=DEVNULL, stderr=DEVNULL) as proc:
        proc.wait()

    return aln_path


def get_all_pham_gene_translations(alchemist):
    """Creates a 2D dictionary that maps phams to dictionaries that map
    unique translations to respective geneids.

    :param alchemist: A connected and fully built AlchemyHandler object
    :type alchemist: AlchemyHandler
    :return: Returns a dictionary mapping phams to translations to geneids
    :rtype: dict{dict}
    """
    engine = alchemist.engine

    # Build phage>>cluster lookup table
    cluster_lookup = dict()
    host_lookup = dict()
    query = "SELECT PhageID, Cluster, Subcluster, HostGenus FROM phage"
    results = engine.execute(query)
    for result in results:
        phageid = result["PhageID"]
        cluster = result["Cluster"]
        subcluster = result["Subcluster"]
        host = result["HostGenus"]
        if cluster is None:
            cluster_lookup[phageid] = "Singleton"
        elif subcluster is None:
            cluster_lookup[phageid] = cluster
        else:
            cluster_lookup[phageid] = subcluster

        host_lookup[phageid] = host

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
            geneid = "".join([geneid, " [product=", str(product), "]"])

        cluster = cluster_lookup[phageid]
        if cluster is not None and cluster > "":
            geneid = "".join([geneid, " [cluster=", str(cluster), "]"])

        host = host_lookup[phageid]
        if host is not None and host > "":
            geneid = "".join([geneid, " [host_genus=", str(host), "]"])

        if locus is not None and locus > "":
            geneid = "".join([geneid, " [locus=", str(locus), "]"])

        pham_translations = phams.get(phamid, dict())
        gene_ids = pham_translations.get(translation, [])
        gene_ids.append(geneid)
        pham_translations[translation] = gene_ids
        phams[phamid] = pham_translations

    engine.dispose()

    return phams


def dump_pham_out_fastas(working_dir, phams_dict, neglect_singles=True,
                         threads=1, verbose=False):
    """Uses multiple threads to write fasta-formatted multiple sequence files
    for all of the phams listed.

    :param working_dir: Path to the directory where the files will be written
    :type working_dir: pathlib.Path
    :param phams_dict: Dictionary that maps phams to translations to geneids
    :type phams_dict: dict{dict}
    :param neglect_singles: Omit single sequence fasta files from the path map
    :type neglect_singles: bool
    :param threads: Number of threads to spawn during export
    :type threads: int
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :return pham_fasta_map: Dictionary that maps phams to their fasta file path
    :rtype pham_fasta_map: dict
    """
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
                            verbose=verbose, lock_timeout=THREAD_LOCK_TIMEOUT,
                            join_timeout=THREAD_JOIN_TIMEOUT)

    return pham_fasta_map


def align_pham_out_fastas(working_dir, pham_fasta_map, threads=1,
                          verbose=False):
    """Uses multiple processes to align fasta-formatted multiple sequence files
    for all of the phams listed.

    :param working_dir: Path to the directory where the files will be written
    :type working_dir: pathlib.Path
    :param phams_dict: Dictionary that maps phams to their fasta file path
    :type phams_dict: dict{Path}
    :param threads: Number of processes/threads to spawn during alignment
    :type threads: int
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :return pham_aln_map: Dictionary that maps phams to their aln file path
    :rtype pham_aln_map: dict
    """
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


def reintroduce_pham_fasta_duplicates(pham_path_map, pham_translations_dict,
                                      threads=1, verbose=False):
    """Uses multiple threads to read and re-write fasta-formatted sequence
    files for all of the phams listed.

    :param pham_path_map: Dictionary that maps phams to their file path
    :type pham_path_map: dict{Path}
    :param phams_translations_dict: Map of phams to translations to geneids
    :type phams_translations_dict: dict{dict}
    :param threads: Number of threads to spawn during reading/writing of files
    :type threads: int
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    """
    work_items = []
    for pham, filepath in pham_path_map.items():
        ts_to_gs = pham_translations_dict.get(pham, {})

        work_items.append((ts_to_gs, filepath))

    multithread.multithread(
                        work_items, threads,
                        fileio.reintroduce_fasta_duplicates, verbose=verbose,
                        lock_timeout=THREAD_LOCK_TIMEOUT,
                        join_timeout=THREAD_JOIN_TIMEOUT)


# Parallelized all-encompassing function that combines the functionality
# of the above functions
def write_phams(fasta_dir, aln_dir, phams_translations_dict, cores=1,
                verbose=False):
    work_items = []
    for pham, pham_translations in phams_translations_dict.items():
        work_items.append((fasta_dir, aln_dir, pham, pham_translations))

    parallelize.parallelize(work_items, cores, write_phams_process,
                            verbose=verbose)


def write_phams_process(fasta_dir, aln_dir, pham, pham_translations):
    fasta_path = fasta_dir.joinpath("".join([str(pham), "_genes.fasta"]))
    aln_path = aln_dir.joinpath("".join([str(pham), "_genes.aln"]))

    gs_to_ts = {}
    for translation, gene_ids in pham_translations.items():
        gs_to_ts[gene_ids[0]] = translation

    fileio.write_fasta(gs_to_ts, fasta_path)

    if len(pham_translations) > 1:
        run_clustalo(fasta_path, aln_path)
        fileio.reintroduce_fasta_duplicates(pham_translations, aln_path)

    fileio.reintroduce_fasta_duplicates(pham_translations, fasta_path)
