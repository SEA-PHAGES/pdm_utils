"""Functions for running rpsblast and parsing the output."""

from Bio.Blast import NCBIXML

from pdm_utils.functions.fasta import write_fasta
from pdm_utils.functions.subprocess import run_command
from pdm_utils.functions.multiprocess import parallelize
from pdm_utils.classes.error import BlastError

BATCH_SIZE = 250


def get_version(program="rpsblast", long=False):
    """Get the version string from BLAST.

    :param program: globally executable rpsblast(+) to get version for
    :type program: str or pathlib.Path
    :param long: return the long version string
    :type long: bool
    :return: version
    """
    command = f"{program} -version"

    try:
        stderr, stdout = run_command(command)
        if stdout:
            raise BlastError(stdout)
        if long:
            version = stderr.split("\n")[1].split("blast ")[-1]
        else:
            version = stderr.split("\n")[0].split(": ")[-1]
    except FileNotFoundError:
        version = None

    return version


def rpsblast(program, query, db, out, evalue, comp_based_stats=1,
             seg="no", outfmt=5):
    """Run rpsblast(+).

    Query the sequences in `infile` against `database` with results
    written to `outfile`.

    :param program: name or path of the rpsblast(+) program to call
    :type program: str or pathlib.Path
    :param query: path to input file containing query sequences
    :type query: str or pathlib.Path
    :param db: path to database containing NCBI conserved domains
    :type db: str or pathlib.Path
    :param out: path to output file to write results to
    :type out: str or pathlib.Path
    :param evalue: significance threshold to retain Cdd hits
    :type evalue: int or float
    :param comp_based_stats: adjust stats for low-complexity regions
    :type comp_based_stats: int
    :param seg: filter query sequence with SEG
    :type seg: str
    :param outfmt: output file options
    :type outfmt: int or str
    :return: out
    :raise: BlastError
    """
    command = f"{program} -query {query} -db {db} -out {out} -outfmt " \
              f"{outfmt} -evalue {evalue} -seg {seg} -comp_based_stats " \
              f"{comp_based_stats}"

    stdout, stderr = run_command(command)
    if stderr:
        # Problem with filepaths or other parameters
        raise BlastError(stderr)

    return out


def blastall(program, sequences, db, tmp_dir, evalue, cpus=1,
             batch_size=BATCH_SIZE):
    """Query each of the given sequences against the NCBI conserved
    domain database using rpsblast.

    Batches are set up to be reasonably load balancing. The more
    sequences are searched, the better the load balancing works.

    :param program: name or path of the rpsblast(+) program to call
    :type program: str or pathlib.Path
    :param sequences: the (geneid, translation) pairs to search
    :type sequences: list[tuple[str, str]]
    :param db: path to database containing NCBI conserve domains
    :type db: str or pathlib.Path
    :param tmp_dir: directory to use for temporary files
    :type tmp_dir: pathlib.Path
    :param evalue: significance threshold to retain Cdd hits
    :type evalue: int or float
    :param batch_size: maximum number of genes to run per batch
    :type batch_size: int
    :param cpus: number of CPU cores to use
    :type cpus: int
    :return: outfiles
    :rtype: list[str or pathlib.Path]
    """
    # Determine appropriate number of partitions for good load balancing
    num_batches = max((len(sequences) // batch_size, 1))
    while len(sequences) / num_batches > batch_size or num_batches % cpus != 0:
        num_batches += 1

    # Sort the sequences by translation length
    sequences = sorted(sequences, key=lambda x: len(x[-1]))

    jobs = list()
    for batch_i in range(num_batches):
        geneids, translations = list(), list()
        for batch_index in range(len(sequences))[batch_i::num_batches]:
            geneid, translation = sequences[batch_index]
            geneids.append(geneid)
            translations.append(translation)

        infile = tmp_dir.joinpath(f"rpsblast_batch_{batch_i}.fasta")
        outfile = tmp_dir.joinpath(f"rpsblast_batch_{batch_i}.xml")

        write_fasta(geneids, translations, infile)

        jobs.append((program, infile, db, outfile, evalue))

    outfiles = parallelize(jobs, cpus, rpsblast, verbose=True)

    return outfiles


def parse_xml(filepath, sequences, evalue):
    """Return conserved domain hits parsed from an XML results file
    output by rpsblast.

    :param filepath: the rpsblast XML output file
    :type filepath: pathlib.Path
    :param sequences: the geneid, translation pairs queried
    :type sequences: list[tuple[str, str]]
    :param evalue: significance threshold to retain HSPs
    :type evalue: float
    :return: domain_data
    """
    sequence_map = dict()
    for geneid, translation in sequences:
        sequence_map[geneid] = translation

    domain_data = dict()
    with open(filepath, "r") as blast_reader:
        for blast_record in NCBIXML.parse(blast_reader):
            query, domains = sequence_map[blast_record.query], list()
            for alignment in blast_record.alignments:
                domain_id, name, description = None, None, None

                definition = alignment.hit_def.replace("\"", "\'")
                def_list = definition.split(",")
                if len(def_list) == 1:
                    description = def_list[0].strip()
                else:
                    domain_id = def_list[0].strip()
                    name = def_list[1].strip()

                    if len(def_list) > 2:
                        description = ",".join(def_list[2:]).strip()

                for hsp in alignment.hsps:
                    if hsp.expect >= evalue:
                        continue

                    domain = {"HitID": alignment.hit_id,
                              "DomainID": domain_id,
                              "Name": name,
                              "Description": description,
                              "Expect": float(hsp.expect),
                              "QueryStart": int(hsp.query_start),
                              "QueryEnd": int(hsp.query_end)}
                    domains.append(domain)

            domain_data[query] = domains

    return domain_data
