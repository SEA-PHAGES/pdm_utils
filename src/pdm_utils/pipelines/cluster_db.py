"""Pipeline for clustering phage genomes in the database."""
import argparse
import csv
import shutil
import string
import sys
import time
from pathlib import Path

from phamclust.matrix import matrix_from_squareform

from pdm_utils.functions import (configfile, fileio, pipelines_basic,
                                 subprocess)
from pdm_utils.pipelines.revise import TICKET_HEADER

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_cluster_db"

TEMP_DIR = "/tmp/pdm_utils_cluster_temp"


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------

def main(unparsed_args_list):
    """Uses parsed args to run the entirety of the file export pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    # Returns after printing appropriate error message from parsing/connecting.
    args = parse_cluster_db(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    values = pipelines_basic.parse_value_input(args.input)
    values = None

    execute_cluster_db(alchemist, folder_path=args.folder_path,
                       folder_name=args.folder_name,
                       values=values, filters=args.filters,
                       verbose=args.verbose, force=args.force,
                       cpus=args.cpus, cluster_size=args.cluster_size,
                       cluster_thres=args.cluster_thres,
                       subcluster_thres=args.subcluster_thres,
                       subcluster=args.subcluster)


def parse_cluster_db(unparsed_args_list):
    """Parses cluster_db arguments and stores them with an argparse object.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    :returns: ArgParse module parsed args.
    """
    DATABASE_HELP = "Name of the MySQL database to export from."

    CONFIG_FILE_HELP = """
        Export option that enables use of a config file for sourcing
        credentials
            Follow selection argument with the path to the config file
            specifying MySQL and NCBI credentials.
        """
    FORCE_HELP = """
        Export option that aggresively creates and overwrites directories.
        """
    VERBOSE_HELP = """
        Export option that enables progress print statements.
        """
    FOLDER_PATH_HELP = """
        Export option to change the path
        of the directory where the exported files are stored.
            Follow selection argument with the path to the
            desired export directory.
        """
    FOLDER_NAME_HELP = """
        Export option to change the name
        of the directory where the exported files are stored.
            Follow selection argument with the desired name.
        """

    IMPORT_FILE_HELP = """
        Selection input option that imports values from a csv file.
            Follow selection argument with path to the
            csv file containing the names of each genome in the first column.
        """
    SINGLE_GENOMES_HELP = """
        Selection input option that imports values from cmd line input.
            Follow selection argument with space separated
            names of genomes in the database.
        """
    WHERE_HELP = """
        Data filtering option that filters data by the inputted expressions.
            Follow selection argument with formatted filter expression:
                {Table}.{Column}={Value}
        """

    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str, help=DATABASE_HELP)

    parser.add_argument("-c", "--config_file",
                        type=pipelines_basic.convert_file_path,
                        help=CONFIG_FILE_HELP)
    parser.add_argument("-m", "--folder_name", type=str,
                        help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path", type=Path,
                        help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)
    parser.add_argument("-f", "--force", action="store_true",
                        help=FORCE_HELP)
    parser.add_argument("--cpus", type=int)

    parser.add_argument("-if", "--import_file",
                        type=pipelines_basic.convert_file_path,
                        help=IMPORT_FILE_HELP, dest="input",
                        default=[])
    parser.add_argument("-in", "--import_names", nargs="*",
                        help=SINGLE_GENOMES_HELP, dest="input")
    parser.add_argument("-w", "--where", nargs="?",
                        help=WHERE_HELP, dest="filters")

    parser.add_argument("-sc", "--subcluster", action="store_true",
                        help="Toggle to perform subclustering")
    parser.add_argument("-k", "--cluster_size", type=int,
                        help="Minimum cluster size to perform subclustering")
    parser.add_argument("-ct", "--cluster_thres", type=float,
                        help="Similarity threshold to use for clustering")
    parser.add_argument("-st", "--subcluster_thres", type=float,
                        help="Similarity threshold to use for subclustering")

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME, folder_path=None,
                        config_file=None, verbose=False, input=[],
                        filters="", cpus=1,
                        cluster_size=6, cluster_thres=0.25,
                        subcluster_thres=0.6)

    parsed_args = parser.parse_args(unparsed_args_list[2:])

    return parsed_args


def execute_cluster_db(alchemist, folder_path=None,
                       folder_name=DEFAULT_FOLDER_NAME, values=None,
                       verbose=False, force=False, filters="", cpus=1,
                       cluster_size=6, cluster_thres=0.25,
                       subcluster=False, subcluster_thres=0.6):
    """Executes the entirety of the phage clustering pipeline.

    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param pipeline: File type that dictates data processing.
    :type pipeline: str
    :param folder_path: Path to a valid dir for new dir creation.
    :type folder_path: Path
    :param folder_name: A name for the export folder.
    :type folder_name: str
    :param force: A boolean to toggle aggresive building of directories.
    :type force: bool
    :param values: List of values to filter database results.
    :type values: list[str]
    :param verbose: A boolean value to toggle progress print statements.
    :param filters: A list of lists with filter values, grouped by ORs.
    :type filters: str
    :param cpus: Number of processes/threads to spawn during the pipeline
    :type cpus: int
    :param k: Number of phage genomes in a cluster required for subclustering.
    :type k: int
    :param c: PEQ threshold between genomes for phage clustering.
    :type c: float
    :param s: PEQ threshold between genomes for phage subclustering.
    :type s: float
    :param subcluster: Toggle pipeline subclustering
    :type subcluster: bool
    """
    # Retrieve information for genomes filtered for
    db_filter = pipelines_basic.build_filter(alchemist, "phage", "",
                                             values=values,
                                             verbose=verbose)
    filtered_hits = set(db_filter.build_values(
                                    where=db_filter.build_where_clauses()))
    if not filtered_hits:
        print("No database entries received from phage.")

    # Retrieve information for all genomes in database
    db_filter = pipelines_basic.build_filter(alchemist, "phage", "",
                                             values=None)
    db_filter.values = db_filter.build_values()

    # Create temporary directories and working directories
    temp_dir = Path(TEMP_DIR)
    if temp_dir.is_dir():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir()
    export_path = pipelines_basic.create_working_path(folder_path, folder_name,
                                                      force=force)
    pipelines_basic.create_working_dir(export_path, force=force)

    if verbose:
        print("Clustering phages in database with phamclust")

    # Dump phage genome information to file and run phamclust
    phamclust_input = dump_genes_and_translations(db_filter, temp_dir)
    phamclust(phamclust_input, export_path, cpus=cpus,
              k=cluster_size, c=cluster_thres, s=subcluster_thres)

    # Find phamclust output directory
    phamclust_dir = None
    for phamclust_dir in export_path.iterdir():
        if not phamclust_dir.is_dir():
            continue

        if phamclust_dir.name == ".DS_Store":
            continue

        if "phamclust" not in phamclust_dir.name:
            continue

        break

    phamclust_dir = export_path # delete me

    if phamclust_dir is None:
        print("ERROR: PhamClust output could not be detected.\n"
              "Please check your PhamClust installation.")
        sys.exit(1)

    if verbose:
        print("Mapping new clusters with the existing cluster schema")
    # Parse phamclust output and assign changed/new cluster assignments
    cluster_scheme, subcluster_scheme = parse_phamclust_cluster_schema(
                                                        phamclust_dir)

    print("cluster_scheme, subcluster_scheme")
    print(cluster_scheme, subcluster_scheme)
    cluster_lookup, seqid_cluster_map, subcluster_lookup = \
        query_cluster_metadata(db_filter)

    old_clusters = list(seqid_cluster_map.keys())
    cluster_redistributions = get_cluster_redistributions(
                                cluster_scheme, cluster_lookup, old_clusters)
    named_cluster_scheme = assign_cluster_names(
                                cluster_scheme, cluster_redistributions,
                                verbose=verbose)
    new_cluster_lookup = dict()
    for cluster, phages in named_cluster_scheme.items():
        for phage in phages:
            new_cluster_lookup[phage] = cluster

    update_dicts = get_cluster_updates(named_cluster_scheme, cluster_lookup)

    if subcluster:
        # Assign changed/new subcluster assignments
        for i, phages in cluster_scheme.items():
            if i is None:
                continue

            cluster = new_cluster_lookup[phages[0]]

            cluster_subcluster_scheme = subcluster_scheme.get(i, dict())
            if len(cluster_subcluster_scheme) <= 1:
                continue

            old_subclusters = set()
            for phage in phages:
                subcluster = subcluster_lookup[phage]

                if subcluster is None:
                    continue

                if cluster == subcluster.rstrip("1234567890"):
                    old_subclusters.add(subcluster)

            cluster_redistributions = get_cluster_redistributions(
                                                cluster_subcluster_scheme,
                                                subcluster_lookup,
                                                old_subclusters)

            named_cluster_subcluster_scheme = assign_cluster_names(
                                                cluster_subcluster_scheme,
                                                cluster_redistributions,
                                                subcluster=cluster,
                                                verbose=verbose)

            update_dicts += get_cluster_updates(
                                        named_cluster_subcluster_scheme,
                                        subcluster_lookup,
                                        field="Subcluster")

    if verbose:
        print("Writing updates to file")
    # Write cluster/subcluster assignments to an update-style csv file
    update_dicts = [update_dict for update_dict in update_dicts
                    if update_dict["key_value"] in filtered_hits]
    print(export_path)
    write_clustering_update_ticket(export_path, update_dicts,
                                   filename="updates.csv")

    # Clean up
    if temp_dir.is_dir():
        shutil.rmtree(temp_dir)


# phamclust helper functions
# =============================================================================
def dump_genes_and_translations(db_filter, output_dir):
    """Dumps PhageID, genes, and their amino acid protein sequences to a
    tab-delimited file for use in PhamClust

    :param db_filter: A connected and fully built Filter object.
    :type db_filter: Filter
    :param output_dir: Path to a dir to dump genes and translations.
    :type output_dir: Path
    :returns: Path to the newly created tab-delimited file
    :rtype: Path
    """

    output_file = output_dir.joinpath("genes.tsv")

    rows = db_filter.select(["gene.PhageID", "gene.PhamID",
                             "gene.Translation"])

    for row in rows:
        row["Translation"] = row["Translation"].decode("utf-8")

    fieldnames = ["PhageID", "PhamID", "Translation"]
    with open(output_file, "w") as filehandle:
        csv_writer = csv.DictWriter(filehandle, fieldnames, delimiter="\t")

        for row in rows:
            csv_writer.writerow(row)

    return output_file


def phamclust(input_file, output_dir, cpus=1, k=6, c=0.25, s=0.6):
    """ Runs phamclust on a tab-delimited file containing phage-encoded
    amino acid sequences and their phamilies.

    :param input_file: Path to an input tab-delimted file.
    :type input_file: Path
    :param output_dir: Path to create an output directory for phamclust.
    :type output_dir: Path
    :param cpus: Number of cores to use during clustering.
    :type cpus: int
    :param k: Number of phage genomes in a cluster required for subclustering.
    :type k: int
    :param c: PEQ threshold between genomes for phage clustering.
    :type c: float
    :param s: PEQ threshold between genomes for phage subclustering.
    :type s: float
    """
    command = (f"phamclust {input_file} {output_dir} -t {cpus} "
               f"-c {c} -s {s} -k {k}")

    subprocess.run_command(command)


def parse_phamclust_cluster_schema(phamclust_dir):
    """ Parses phamclust output for phage genome clusters and subclusters.
    :param phamclust_dir: Path to an output directory produced by phamclust.
    :type phamclust_dir: Path
    :returns: Dicts containing information mapping clusters to phage genomes.
    :rtype: (dict, dict).
    """
    # Initiate schema dictionaries
    print("phamclust_dir:", phamclust_dir)

    # Cluster schema is a dict where keys are valid cluster indicies and
    # values are lists of phage genome names
    cluster_schema = dict()

    # Subcluster schema is a dict where keys are valid cluster indicies,
    # values are also dicts containing keys with valid subcluster indicies
    # and values are lists of phage genome names
    subcluster_schema = dict()

    for cluster_dir in phamclust_dir.iterdir():
        print("cluster_dir:", cluster_dir)

        # Skip Mac associated files and other nonsense
        if not cluster_dir.is_dir():
            continue

        if cluster_dir.name == ".DS_Store":
            continue

        # Parse phage genome names of singleons (unclustered phage genomes)
        # by iterating through unclustered phage genome files
        if "singleton" in cluster_dir.name:
            genomes_dir = cluster_dir.joinpath("genomes")
            genomes = list()

            for genome_file in genomes_dir.iterdir():
                if genome_file.suffix != ".faa":
                    continue

                genomes.append(genome_file.stem)

            cluster_schema[None] = genomes

        # Guess that all phage genome clusters will have an output folder
        # labelled as cluster_*
        if "cluster" not in cluster_dir.name:
            continue

        # Parse cluster index from cluster_*
        cluster_num = int(cluster_dir.name.split("_")[1])

        # Iterate through .tsv files containing distance matricies for the
        # total cluster and for subclusters
        cluster_subcluster_schema = dict()
        for dist_file in cluster_dir.iterdir():
            if dist_file.suffix != ".tsv":
                continue

            sym_matrix = matrix_from_squareform(dist_file)
            genomes = sym_matrix.nodes

            # Total cluster distance matrix is usually peq_similarity.tsv
            if "peq" in dist_file.name:
                cluster_genomes = genomes
                continue

            # Subcluster-specific distance matricies are assumed to be
            # subcluster_*_similarity.tsv
            if "subcluster" not in dist_file.name:
                continue

            subcluster_num = int(dist_file.name.split("_")[1])
            cluster_subcluster_schema[subcluster_num] = genomes

        cluster_schema[cluster_num] = cluster_genomes
        subcluster_schema[cluster_num] = cluster_subcluster_schema

    return cluster_schema, subcluster_schema


# Cluster naming helper functions
# =============================================================================
def assign_cluster_names(cluster_scheme, cluster_redistributions,
                         verbose=False, cluster_prefix=None, subcluster=None):
    """ Assigns cluster/subcluster names given the information about the
    current schema. Cluster names are generated where possible, with cluster
    names receiving the next available character(s) and subclusters receiving
    concatenations of the cluster and their subcluster index.

    :param cluster_scheme: A dict mapping cluster indicies to phage genomes.
    :type cluster_scheme: dict
    :param cluster_redistributions: Phage genomes which have changed clusters
    :type cluster_redistributions: dict
    :param verbose: Toggle verbosity.
    :type verbose: bool
    :param cluster_prefix: Prefix to append to given clusters.
    :type cluster_prefix: str
    :param subcluster: Specify cluster name for subcluster assignment
    :type subcluster: str
    :returns: A dict mapping cluster names to phage genome names.
    :rtype: dict
    """
    named_scheme = dict()
    named_scheme[None] = cluster_scheme.pop(None, list())

    old_clusters = list(cluster_redistributions.keys())

    remaining = list([x for x in cluster_scheme.keys() if x is not None])
    assigned = list()

    for old_cluster, cluster_redistribution in cluster_redistributions.items():
        if old_cluster is None:
            continue

        old_cluster_inheretors = sorted(cluster_redistribution.items(),
                                        key=lambda x: len(x[1]), reverse=True)

        for old_cluster_inheretor in old_cluster_inheretors:
            num_cluster = old_cluster_inheretor[0]
            if num_cluster is None:
                continue

            if num_cluster in assigned:
                continue

            named_scheme[old_cluster] = cluster_scheme[num_cluster]
            remaining.remove(num_cluster)
            assigned.append(num_cluster)
            break

    for num_cluster in remaining:
        if subcluster is not None:
            if len(cluster_scheme) <= 1:
                named_scheme[None] = cluster_scheme[num_cluster]
                continue

            new_subcluster = gen_new_subcluster(subcluster, old_clusters)
            named_scheme[new_subcluster] = cluster_scheme[num_cluster]

            if verbose:
                print(f"......Created new subcluster '{new_subcluster}'...")

            old_clusters.append(new_subcluster)
            assigned.append(num_cluster)
            continue

        new_cluster = gen_new_cluster(old_clusters,
                                      cluster_prefix=cluster_prefix)
        named_scheme[new_cluster] = cluster_scheme[num_cluster]

        old_clusters.append(new_cluster)
        assigned.append(num_cluster)

        if verbose:
            print(f"...Created new cluster '{new_cluster}'...")

    return named_scheme


def gen_new_cluster(old_clusters, cluster_prefix=None):
    """Generate a new cluster name with capital alphabetic characters.
    Clusters which have exceeded the traditional english alphabet will
    have additional significant characters added to the cluster i.e.
    after cluster Z is cluster AA.

    :param old_clusters: Used cluster names
    :type old_clusters: list
    :param cluster_prefix: Prefix to append to the front of a cluster name.
    :type cluster_prefix: str
    :returns: Newly generated cluster
    :rtype: str
    """
    if cluster_prefix is None:
        cluster_prefix = ""
    alphabet = string.ascii_uppercase

    chars = 1
    counter = []
    while True:
        filled = True
        if chars > len(counter):
            for _ in range(chars - len(counter)):
                counter.append(0)

        letters = [cluster_prefix]
        for char_num in counter:
            letters.append(alphabet[char_num])

        new_cluster = "".join(letters)
        if new_cluster not in old_clusters:
            return new_cluster

        for i in range(len(counter)):
            if counter[len(counter)-1-i] < (len(alphabet) - 1):
                counter[len(counter)-1-i] += 1
                filled = False
                break

            counter[len(counter)-1-i] = 0

        if filled:
            chars += 1
            if chars >= 5:
                break


def gen_new_subcluster(cluster, old_subclusters):
    for i in range(1, 1000):
        subcluster = "".join([cluster, str(i)])
        if subcluster not in old_subclusters:
            return subcluster


# Reclustering analysis helper functions
# =============================================================================

def query_cluster_metadata(db_filter):
    """ Retrieve cluster/subcluster metadata.

    :param db_filter: A connected Filter object:
    :type db_filter: pdm_utils.classes.Filter
    :returns: Lookup dictionaries for clusters and subclusters
    :rtype: (dict, dict, dict)
    """
    cluster_data = db_filter.retrieve(["phage.Cluster",
                                       "phage.Subcluster"])
    cluster_lookup = {}
    subcluster_lookup = {}
    for phage_id, data_dict in cluster_data.items():
        cluster_lookup[phage_id] = data_dict["Cluster"][0]
        subcluster_lookup[phage_id] = data_dict["Subcluster"][0]

    seqid_cluster_map = db_filter.group("phage.Cluster")

    return (cluster_lookup, seqid_cluster_map,
            subcluster_lookup)


def get_cluster_redistributions(cluster_scheme, cluster_lookup, old_clusters):
    """ Maps and returns phages which have been reassigned clusters
    (when phage clusters merge or are split.)

    :param cluster_scheme: A dict mapping cluster names to phage genome names.
    :type cluster_scheme: dict
    :param cluster_lookup: A dict mapping phage genome names to old clusters.
    :type cluster_lookup: dict
    :param old_clusters: A set containing all cluster names recorded.
    :type old_clusters: set
    :returns: A map of cluster reassignments.
    :rtype: dict
    """
    cluster_redistributions = dict()

    for cluster_num, cluster_members in cluster_scheme.items():
        for member in cluster_members:
            old_cluster = cluster_lookup.get(member)

            cluster_redistribution = cluster_redistributions.get(
                                                        old_cluster, dict())

            distributed_members = cluster_redistribution.get(
                                                        cluster_num, list())
            distributed_members.append(member)

            cluster_redistribution[cluster_num] = distributed_members
            cluster_redistributions[old_cluster] = cluster_redistribution

    cluster_redistribution_weights = list()
    for old_cluster, cluster_redistribution in cluster_redistributions.items():
        weight = 0
        for new_cluster, cluster_members in cluster_redistribution.items():
            weight += len(cluster_members)

        cluster_redistribution_weights.append((old_cluster, weight))

    cluster_redistribution_weights.sort(key=lambda x: x[1], reverse=True)

    cluster_redistributions = {
                cluster: cluster_redistributions[cluster]
                for cluster, weight in cluster_redistribution_weights}

    return cluster_redistributions


def diff_cluster_schemes(cluster_scheme, cluster_lookup):
    """ Maps and returns all changes created to the clustering schema.

    :param cluster_scheme: A dict mapping cluster names to phage genome names.
    :type cluster_scheme: dict
    :param cluster_lookup: A dict mapping phage genome names to old clusters.
    :type cluster_lookup: dict
    :returns: A map of clustering schema changes.
    :rtype: dict
    """
    scheme_alterations = dict()

    for new_cluster, cluster_members in cluster_scheme.items():
        for member in cluster_members:
            old_cluster = cluster_lookup[member]
            if new_cluster != old_cluster:
                altered_data = scheme_alterations.get(new_cluster, list())
                data = dict()
                data["id"] = member
                data["old_cluster"] = old_cluster
                data["new_cluster"] = new_cluster

                altered_data.append(data)
                scheme_alterations[new_cluster] = altered_data

    return scheme_alterations


def get_cluster_updates(cluster_scheme, cluster_lookup, field="Cluster"):
    """ Convert changes to a clustering scheme into dictionaries
    for database updates.

    :param cluster_scheme: A dict mapping cluster names to phage genome names.
    :type cluster_scheme: dict
    :param cluster_lookup: A dict mapping phage genome names to old clusters.
    :type cluster_lookup: dict
    :param field: The field to update.
    :returns: A list of update tickets
    :rtype: list
    """
    print(cluster_scheme, cluster_lookup)

    scheme_alterations = diff_cluster_schemes(cluster_scheme, cluster_lookup)
    print(scheme_alterations)

    update_dicts = []

    for cluster, diff_data in scheme_alterations.items():
        for data_dict in diff_data:
            update_dict = {}

            if cluster is None:
                cluster = "NULL"

            update_data = ("phage", field, cluster, "PhageID", data_dict["id"])
            for i in range(len(TICKET_HEADER)):
                update_dict[TICKET_HEADER[i]] = update_data[i]

            update_dicts.append(update_dict)

    return update_dicts


def write_clustering_update_ticket(working_dir, update_dicts,
                                   filename=None):
    """ Write update tickets to file.

    :param working_dir: Directory to write file in
    :type working_dir: Path
    :param update_dicts: List of update tickets.
    :type update_dicts: list
    :param filename: Name of the file to be written/
    :type filename: str
    """
    if filename is None:
        filename = working_dir.with_suffix(".csv").name

    if not update_dicts:
        print("no update_dicts")
        return False

    filepath = working_dir.joinpath(filename)
    fileio.export_data_dict(update_dicts, filepath, TICKET_HEADER,
                            include_headers=True)

    return True
