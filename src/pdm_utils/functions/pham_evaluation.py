import shlex
from subprocess import Popen, PIPE

from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes.clustal import *
from pdm_utils.functions import mysqldb_basic


def get_pham_genes(engine, phamid):
    """
    Queries the database for the geneids and translations found in the
    indicated pham. Returns a dictionary mapping the pham's geneids to
    their associated translations. All geneid:translation pairs will be
    represented (i.e. no redundant gene filtering is done...).
    :param engine: the Engine allowing access to the database
    :type engine: sqlalchemy Engine
    :param phamid: the pham whose genes are to be returned
    :type phamid: int
    :return: pham_genes
    :rtype: dict
    """
    pham_genes = dict()

    # Query will return the pham's GeneIDs and Translations grouped by genes
    # that share the same sequence
    query = f"SELECT GeneID, Translation FROM gene WHERE PhamID = {phamid} " \
            f"ORDER BY Translation, GeneID ASC"
    query_results = mysqldb_basic.query_dict_list(engine, query)

    for dictionary in query_results:
        geneid = dictionary["GeneID"]
        translation = dictionary["Translation"].decode("utf-8")
        pham_genes[geneid] = translation

    return pham_genes


def write_genes_to_fasta(genes, infile):
    """
    Writes the input genes to the indicated file in FASTA multiple
    sequence format (unaligned).
    :param genes: the genes to be written to file
    :type genes: dict
    :param infile: the name of the file to write the genes to
    :type infile: str
    :return:
    """
    with open(infile, "w") as fh:
        for geneid, translation in genes.items():
            fh.write(f">{geneid}\n{translation}\n")


def clustalo(fasta_file, aln_out=None, mat_out=None, threads=1, verbose=0):
    """
    Runs Clustal Omega to generate a multiple sequence alignment (MSA)
    and percent identity matrix (PIM) for the indicated file. Infile is
    expected to be in FASTA multiple sequence format. MSA will be in
    Clustal format.
    :param fasta_file: FASTA file containing sequences to be aligned
    :type fasta_file: str
    :param aln_out: the multiple sequence alignment (MSA) output file
    :type aln_out: str
    :param mat_out: the percent identity matrix (PIM) output file
    :type mat_out: str
    :param threads: number of threads to use for alignment
    :type threads: int
    :param verbose: verbosity level (0-2)
    :type verbose: int
    :return: [aln_out, mat_out]
    :rtype: list
    """
    # Create sensible output filenames if none are provided
    basename = fasta_file.split(".")[0]
    if aln_out is None:
        aln_out = f"{basename}.aln"
    if mat_out is None:
        mat_out = f"{basename}.mat"

    # Make sure verbose is in range(0,3)
    if verbose <= 0:
        verbose = 0
    elif verbose > 2:
        verbose = 2

    # Build Clustal Omega command that will produce a clustal-formatted
    # alignment output file and percent identity matrix
    command = f"clustalo -i {fasta_file} -o {aln_out} --distmat-out=" \
              f"{mat_out} --outfmt=clustal --full --percent-id --force " \
              f"--output-order=tree-order --threads={threads}"
    for x in range(verbose):
        command += " -v"                # Add verbosity to command
    command = shlex.split(command)      # Convert command to arg list

    # Run the Clustal Omega command as a subprocess
    with Popen(args=command, stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    # Return output files so user can find them easily, whether they specified
    # the filenames or not
    return [aln_out, mat_out]


#################
# Example Usage #
#################
def main():
    # Set up MySQL handler
    alchemist = AlchemyHandler(database="Actinobacteriophage")
    alchemist.connect(pipeline=True)
    engine = alchemist.engine

    # Choose a pham, set up infile name, and construct FASTA
    pham = 49991
    infile = "/tmp/pham_49991.fasta"
    genes = get_pham_genes(engine, pham)
    write_genes_to_fasta(genes, infile)

    # Run clustal and expand the return value to output filenames
    out_aln, out_mat = clustalo(infile)

    # Parse the identity matrix, get the centroid name and average
    # distance from centroid to each gene in the pham
    mat = PercentIdentityMatrix(out_mat)
    mat.parse_matrix()
    centroid = mat.get_centroid()
    average_dist = round(100 - mat.get_average_identity(centroid), 3)

    # Parse the multiple sequence alignment and get centroid sequence w/o gaps
    aln = MultipleSequenceAlignment(out_aln)
    aln.parse_alignment()
    centroid_seq = aln.get_sequence(centroid, gaps=False)

    print("Summary")
    print("=======================================")
    print(f"Centroid gene: {centroid}")
    print(f"Centroid sequence: {centroid_seq[:16]}...")
    print(f"Average distance from centroid: {average_dist}%")


if __name__ == "__main__":
    main()
