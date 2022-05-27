"""Functions that are used in the phameration pipeline"""

import random
import colorsys

from pdm_utils.functions import mysqldb
from pdm_utils.functions import mysqldb_basic


# DATABASE FUNCTIONS
def get_pham_geneids(engine):
    """
    Queries the database for those genes that are already phamerated.
    :param engine: the Engine allowing access to the database
    :return: pham_geneids
    """
    pham_geneids = dict()

    geneid_query = "SELECT GeneID, PhamID FROM gene WHERE PhamID IS NOT NULL"
    geneid_results = mysqldb_basic.query_dict_list(engine, geneid_query)

    for dictionary in geneid_results:
        pham_id = dictionary["PhamID"]
        geneid = dictionary["GeneID"]

        if pham_id in pham_geneids.keys():
            pham_geneids[pham_id] = pham_geneids[pham_id] | {geneid}
        else:
            pham_geneids[pham_id] = {geneid}

    return pham_geneids


def get_pham_colors(engine):
    """
    Queries the database for the colors of existing phams
    :param engine: the Engine allowing access to the database
    :return: pham_colors
    """
    pham_colors = dict()

    color_query = "SELECT PhamID, Color FROM pham"
    color_results = mysqldb_basic.query_dict_list(engine, color_query)

    for dictionary in color_results:
        pham_id = dictionary["PhamID"]
        color = dictionary["Color"]

        pham_colors[pham_id] = color.upper()

    return pham_colors


def get_new_geneids(engine):
    """
    Queries the database for those genes that are not yet phamerated.
    :param engine: the Engine allowing access to the database
    :return: new_geneids
    """
    new_geneids = list()

    gene_query = "SELECT GeneID FROM gene WHERE PhamID IS NULL"
    gene_results = mysqldb_basic.query_dict_list(engine, gene_query)

    # At scale, much cheaper to convert list to set than to build the
    # set one gene at a time
    for dictionary in gene_results:
        geneid = dictionary["GeneID"]
        new_geneids.append(geneid)

    return set(new_geneids)


def get_geneids_and_translations(engine):
    """
    Constructs a dictionary mapping all geneids to their translations.
    :param engine: the Engine allowing access to the database
    :return: gs_to_ts
    """
    gs_to_ts = dict()

    query = ("SELECT GeneID, CONVERT(Translation USING utf8) as Translation "
             "FROM gene")
    results = mysqldb_basic.query_dict_list(engine, query)

    for dictionary in results:
        geneid = dictionary["GeneID"]
        translation = dictionary["Translation"]
        gs_to_ts[geneid] = translation

    return gs_to_ts


def get_translation_groups(engine):
    """
    Constructs a dictionary mapping all unique translations to their
    groups of geneids that share them
    :param engine: the Engine allowing access to the database
    :return: ts_to_gs
    """
    ts_to_gs = dict()

    query = ("SELECT GeneID, CONVERT(Translation USING utf8) as Translation "
             "FROM gene")
    results = mysqldb_basic.query_dict_list(engine, query)

    for dictionary in results:
        geneid = dictionary["GeneID"]
        trans = dictionary["Translation"]
        geneids = ts_to_gs.get(trans, [])
        geneids.append(geneid)
        ts_to_gs[trans] = geneids

    return ts_to_gs


def update_pham_table(colors, engine):
    """
    Populates the pham table with the new PhamIDs and their colors.
    :param colors: new pham color data
    :type colors: dict
    :param engine: sqlalchemy Engine allowing access to the database
    :return:
    """
    # First command needs to clear the pham table
    commands = ["DELETE FROM pham"]

    # Then we need to issue insert commands for each pham
    for key in colors.keys():
        commands.append(f"INSERT INTO pham (PhamID, Color) VALUES ({key}, "
                        f"'{colors[key]}')")

    mysqldb.execute_transaction(engine, commands)


def update_gene_table(phams, engine):
    """
    Updates the gene table with new pham data
    :param phams: new pham gene data
    :type phams: dict
    :param engine: sqlalchemy Engine allowing access to the database
    :return:
    """
    commands = []

    # We need to issue an update command for each gene in each pham
    for key in phams.keys():
        for gene in phams[key]:
            commands.append(f"UPDATE gene SET PhamID = {key} WHERE GeneID = '{gene}'")

    mysqldb.execute_transaction(engine, commands)


def fix_white_phams(engine):
    """
    Find any phams with 2+ members which are colored as though they are
    orphams (#FFFFFF in pham.Color).
    :param engine: sqlalchemy Engine allowing access to the database
    :return:
    """
    query = "SELECT c.PhamID FROM (SELECT g.PhamID, COUNT(GeneID) AS count, "\
            "p.Color FROM gene AS g INNER JOIN pham AS p ON g.PhamID " \
            "= p.PhamID GROUP BY PhamID) AS c WHERE Color = '#FFFFFF' "\
            "AND count > 1"

    results = mysqldb_basic.query_dict_list(engine, query)
    print(f"Found {len(results)} white phams...")

    commands = []
    for dictionary in results:
        pham_id = dictionary["PhamID"]
        h = s = v = 0
        while h <= 0:
            h = random.random()
        while s < 0.5:
            s = random.random()
        while v < 0.8:
            v = random.random()
        rgb = colorsys.hsv_to_rgb(h, s, v)
        rgb = (rgb[0] * 255, rgb[1] * 255, rgb[2] * 255)
        hexrgb = "#{:02x}{:02x}{:02x}".format(int(rgb[0]), int(rgb[1]),
                                              int(rgb[2]))
        new_color = hexrgb.upper()
        commands.append(f"UPDATE pham SET Color = '{new_color}' WHERE "
                        f"PhamID = '{pham_id}'")

    mysqldb.execute_transaction(engine, commands)


def fix_colored_orphams(engine):
    """
    Find any single-member phams which are colored as though they are
    multi-member phams (not #FFFFFF in pham.Color).
    :param engine: sqlalchemy Engine allowing access to the database
    :return:
    """
    query = "SELECT * FROM (SELECT g.PhamID, COUNT(GeneID) AS count, " \
            "p.Color FROM gene AS g INNER JOIN pham AS p ON g.PhamID " \
            "=p.PhamID GROUP BY PhamID) AS c WHERE Color != '#FFFFFF' " \
            "AND count = 1"

    results = mysqldb_basic.query_dict_list(engine, query)
    print(f"Found {len(results)} non-white orphams...")

    commands = []
    for dictionary in results:
        pham_id = dictionary["PhamID"]
        count = dictionary["count"]
        color = dictionary["Color"]
        new_color = "#FFFFFF"
        commands.append(f"UPDATE pham SET Color = '{new_color}' WHERE "
                        f"PhamID = '{pham_id}'")

    mysqldb.execute_transaction(engine, commands)


def preserve_phams(old_phams, new_phams, old_colors, new_genes):
    """
    Attempts to keep pham numbers consistent from one round of pham
    building to the next
    :param old_phams: the dictionary that maps old phams to their genes
    :param new_phams: the dictionary that maps new phams to their genes
    :param old_colors: the dictionary that maps old phams to colors
    :param new_genes: the set of previously un-phamerated genes
    :return:
    """
    final_phams = dict()
    final_colors = dict()

    outcount = 0
    total = len(old_phams)
    new_phams_copy = new_phams.copy()

    # Iterate through old and new phams
    for old_key in old_phams.keys():
        outcount += 1
        # print("Pham Name Conservation: {} / {}".format(outcount, total))

        old_pham = old_phams[old_key]

        if old_key in final_phams.keys():
            continue

        for new_key in new_phams_copy.keys():
            new_pham = new_phams_copy[new_key]

            if old_pham & new_pham == set():
                continue

            # Case 1 + 5 (Identity and Subtraction)
            if old_pham == new_pham:
                final_phams[old_key] = new_pham
                final_colors[old_key] = old_colors[old_key]
                new_phams_copy.pop(new_key)
                break

            # Case 2 and 4 (Addition and Join) - PHAM GREW
            elif new_pham - old_pham != set():

                # Case 2 and 4 (Addition and Join)
                if new_pham & new_genes != set():

                    # Case 4 - Join with new gene
                    if (new_pham - (new_pham & new_genes)) - old_pham != set():
                        break

                    # Case 2 - Addition with new gene
                    final_phams[old_key] = new_pham
                    final_colors[old_key] = old_colors[old_key]
                    new_phams_copy.pop(new_key)
                    break

                # Case 4 - Join without new gene
                else:
                    break

            # Case 3 - split - PHAM SHRANK, BUT NOT BY REMOVAL
            elif old_pham - new_pham != set():
                break

    final_phams[0] = "placeholder"
    highest_pham = max(map(int, final_phams.keys())) + 1
    final_phams.pop(0)

    # Reassign data for split or joined phams
    for key in new_phams_copy.keys():
        new_key = highest_pham
        highest_pham += 1

        final_phams[new_key] = new_phams_copy[key]

        # Multi-member pham gets non-white color
        if len(new_phams_copy[key]) > 1:
            h = s = v = 0
            while h <= 0:
                h = random.random()
            while s < 0.5:
                s = random.random()
            while v < 0.8:
                v = random.random()
            rgb = colorsys.hsv_to_rgb(h, s, v)
            rgb = (rgb[0] * 255, rgb[1] * 255, rgb[2] * 255)
            hexrgb = "#{:02x}{:02x}{:02x}".format(int(rgb[0]), int(rgb[1]),
                                                  int(rgb[2]))
            final_colors[new_key] = hexrgb.upper()
        # Orpham gets white
        else:
            final_colors[new_key] = '#FFFFFF'

    return final_phams, final_colors
