"""Mock genome data for building mock databases and mock flat files."""
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqFeature import ExactPosition, BeforePosition, AfterPosition, Reference

unittest_file = Path(__file__)
test_dir = unittest_file.parent
test_file_dir = Path(test_dir, "test_files")
base_flat_file = Path("test_flat_file_10.gb") # Alice genome
base_flat_file_path = Path(test_file_dir, base_flat_file)

# Format of the date the script imports into the database.
current_date = datetime.today().replace(hour=0, minute=0,
                                        second=0, microsecond=0)


###Generic BioPython object instantiation functions

def get_seq(filepath):
    """Get genome sequence from a flat file so that it can be added
    to a MySQL database for testing."""
    seqrecord = SeqIO.read(filepath, "genbank")
    return seqrecord.seq

def create_1_part_seqfeature(start=0, stop=0, strand=1, type="",
                             fuzzy="neither", qualifiers=None):
    """Constructs simple BioPython SeqFeature.

    Start = int
    Stop = int
    Strand = int (-1, 1)
    Type = 'CDS', 'Source', 'tRNA', etc.
    Fuzzy = 'start', 'stop', 'both', or 'neither'
    Qualifiers = dictionary of feature descriptions."""
    if fuzzy == "start":
        seq_ftr = SeqFeature(
                    FeatureLocation(
                        BeforePosition(start),
                        ExactPosition(stop),
                        strand=strand),
                    type=type,
                    qualifiers=qualifiers)
    elif fuzzy == "stop":
        seq_ftr = SeqFeature(
                    FeatureLocation(
                        ExactPosition(start),
                        AfterPosition(stop),
                        strand=strand),
                    type=type,
                    qualifiers=qualifiers)
    elif fuzzy == "both":
        seq_ftr = SeqFeature(
                    FeatureLocation(
                        BeforePosition(start),
                        AfterPosition(stop),
                        strand=strand),
                    type=type,
                    qualifiers=qualifiers)
    else:
        seq_ftr = SeqFeature(
                    FeatureLocation(
                        ExactPosition(start),
                        ExactPosition(stop),
                        strand=strand),
                    type=type,
                    qualifiers=qualifiers)
    return seq_ftr

def create_2_part_seqfeature(start1=0, stop1=0, strand1=1,
                             start2=0, stop2=0, strand2=1, type="",
                             fuzzy="neither", qualifiers=None):
    """Constructs simple BioPython SeqFeature.

    Start1, Start2 = int
    Stop1, Stop2 = int
    Strand1, Strand2 = int (-1, 1)
    Type = 'CDS', 'Source', 'tRNA', etc.
    Fuzzy = 'start', or 'neither'
    Qualifiers = dictionary of feature descriptions."""
    # This function could be improved if needed, but right now only
    # can toggle one of the coordinate fuzziness.
    if fuzzy == "start":
        seq_ftr = SeqFeature(
                    CompoundLocation(
                        [FeatureLocation(
                            BeforePosition(start1),
                            ExactPosition(stop1),
                            strand=strand1),
                        FeatureLocation(
                            ExactPosition(start2),
                            ExactPosition(stop2),
                            strand=strand2)],
                        "join"),
                    type=type,
                    location_operator="join",
                    qualifiers=qualifiers)
    else:
        seq_ftr = SeqFeature(
                    CompoundLocation(
                        [FeatureLocation(
                            ExactPosition(start1),
                            ExactPosition(stop1),
                            strand=strand1),
                        FeatureLocation(
                            ExactPosition(start2),
                            ExactPosition(stop2),
                            strand=strand2)],
                        "join"),
                    type=type,
                    location_operator="join",
                    qualifiers=qualifiers)
    return seq_ftr

def create_3_part_seqfeature(start1=0, stop1=0, strand1=1,
                             start2=0, stop2=0, strand2=1,
                             start3=0, stop3=0, strand3=1, type="",
                             qualifiers=None):
    """Constructs simple BioPython SeqFeature.

    Start1, Start2, Start3 = int
    Stop1, Stop2, Stop3 = int
    Strand1, Strand2, Strand3 = int (-1, 1)
    Type = 'CDS', 'Source', 'tRNA', etc.
    Qualifiers = dictionary of feature descriptions."""
    seq_ftr = SeqFeature(
                CompoundLocation(
                    [FeatureLocation(
                        ExactPosition(start1),
                        ExactPosition(stop1),
                        strand=strand1),
                    FeatureLocation(
                        ExactPosition(start2),
                        ExactPosition(stop2),
                        strand=strand2),
                    FeatureLocation(
                        ExactPosition(start3),
                        ExactPosition(stop3),
                        strand=strand3)],
                    "join"),
                type=type,
                location_operator="join",
                qualifiers=qualifiers)
    return seq_ftr

def create_reference(author_string=None):
    """Returns mock Reference data."""
    reference = Reference()
    reference.authors = author_string
    return reference




###Alice genome

def get_alice_genome_draft_data_in_db(seq=True):
    """Returns mock data from phage table for 'draft' Alice genome.

    Seq = True/False, enables retrieving entire genome sequence if needed."""
    if seq:
        genome_seq = get_seq(base_flat_file_path)
    else:
        genome_seq = ""
    data_dict = {
        "PhageID": "Alice",
        "Accession": "",
        "Name": "Alice_Draft",
        "HostGenus": "Mycobacterium",
        "Sequence": genome_seq,
        "Length": 153401,
        "GC": 64.6808,
        "Status": "draft",
        "DateLastModified": current_date,
        "RetrieveRecord": 1,
        "AnnotationAuthor": 1,
        "Cluster": "C",
        "Subcluster": "C1",
        "Notes": "NULL"
        }
    return data_dict

def get_alice_genome_final_data_in_db():
    """Returns mock data from phage table for 'final' Alice genome."""
    data_dict = get_alice_genome_draft_data_in_db()
    data_dict["Name"] = "Alice"
    data_dict["Accession"] = "JF704092"
    data_dict["Status"] = "final"
    data_dict["DateLastModified"] = current_date
    return data_dict


# Alice CDS 139
# Example 1-part feature on bottom strand.

alice_cds_139_translation = (
    "MKKLVTGGLVAAGITLGLISAPSASAEYMTICPSQVSAVVTANTSCGFADNV"
    "FRGFYRQSGWDPLAYSPATGKVYRMHCAPATTTNWGEAKRCWGVGYGGDLLVVYID"
    )

def get_alice_cds_139_draft_data_in_db():
    """Returns mock data from gene table for 'draft' Alice CDS 139."""
    dict = {
        "GeneID": "Alice_CDS_2",
        "PhageID": "Alice",
        "Start": 88120,
        "Stop": 88447,
        "Parts": 1,
        "Length": (len(alice_cds_139_translation) * 3) + 3,
        "Name": "139",
        "Translation": alice_cds_139_translation,
        "Orientation": "R",
        "Notes": "",
        "LocusTag": None,
        "DomainStatus": 0,
        "PhamID": "NULL"
        }
    return dict

def get_alice_cds_139_qualifier_dict():
    """Returns data to construct Alice CDS 139 flat file feature."""
    dict = OrderedDict(
            [("gene", ["139"]),
             ("locus_tag", ["ALICE_139"]),
             ("note", ["gp139"]),
             ("codon_start", ["1"]),
             ("transl_table", ["11"]),
             ("product", ["hypothetical protein"]),
             ("protein_id", ["AEJ94477.1"]),
             ("translation", [])])
    return dict

def get_alice_cds_139_seqfeature():
    """Constructs Alice CDS 139 flat file feature."""
    seq_ftr = create_1_part_seqfeature(88120, 88447, -1, "CDS")
    return seq_ftr




# Alice CDS 193
# Example 1-part feature on top strand.

alice_cds_193_translation = (
    "MGNNRPTTRTLPTGEKATHHPDGRLVLKPKNSLADALSGQLDAQQEASQALTEALV"
    "QTAVTKREAQADGNPVPEDKRVF"
    )

def get_alice_cds_193_draft_data_in_db():
    """Returns mock data from gene table for 'draft' Alice CDS 193."""
    dict = {
        "GeneID": "Alice_CDS_3",
        "PhageID": "Alice",
        "Start": 110297,
        "Stop": 110537,
        "Parts": 1,
        "Length": (len(alice_cds_193_translation) * 3) + 3,
        "Name": "193",
        "Translation": alice_cds_193_translation,
        "Orientation": "F",
        "Notes": "",
        "LocusTag": None,
        "DomainStatus": 0,
        "PhamID": "NULL"
        }
    return dict

def get_alice_cds_193_seqfeature():
    """Constructs Alice CDS 193 flat file feature."""
    seq_ftr = create_1_part_seqfeature(110297, 110537, 1, "CDS")
    return seq_ftr

def get_alice_cds_193_qualifier_dict():
    """Returns data to construct Alice CDS 193 flat file feature."""
    dict = OrderedDict(
            [("gene", ["193"]),
             ("locus_tag", ["ALICE_193"]),
             ("note", ["gp193"]),
             ("codon_start", ["1"]),
             ("transl_table", ["11"]),
             ("product", ["hypothetical protein"]),
             ("protein_id", ["AEJ94419.1"]),
             ("translation", [])])
    return dict




# Alice CDS 252
# Example 2-part compound feature on top strand that wraps around genome termini.

alice_cds_252_translation = (
    "MFDNHPLMSTLITSVRDLNRLANGVIIRNRCEKCDTSAGPDQGV"
    "HFVKTPLGWLYSDPKGKPTTWELFPSDAIHLPVQVIHRVSESEDVSEITENSPHTRKD"
    "DSEATEGAAPSFEPLYPSSHKTERFTVQDTPDGLGAYVFDFGGDAFGAQTCADALAKV"
    "TGKTWYVMHKTVVENTGIYSSMVRHEEESSS")

def get_alice_cds_252_draft_data_in_db():
    """Returns mock data from gene table for 'draft' Alice CDS 252."""
    dict = {
        "GeneID": "Alice_CDS_4",
        "PhageID": "Alice",
        "Start": 152829,
        "Stop": 4,
        "Parts": 2,
        "Length": (len(alice_cds_252_translation) * 3) + 3,
        "Name": "252",
        "Translation": alice_cds_252_translation,
        "Orientation": "F",
        "Notes": "",
        "LocusTag": None,
        "DomainStatus": 0,
        "PhamID": "NULL"
        }
    return dict

def get_alice_cds_252_qualifier_dict():
    """Returns data to construct Alice CDS 252 flat file feature."""
    dict = OrderedDict(
            [("gene", ["252"]),
             ("locus_tag", ["ALICE_252"]),
             ("note", ["gp252"]),
             ("codon_start", ["1"]),
             ("transl_table", ["11"]),
             ("product", ["hypothetical protein"]),
             ("protein_id", ["AEJ94474.1"]),
             ("translation", [])])
    return dict

def get_alice_cds_252_seqfeature():
    """Constructs Alice CDS 252 flat file feature."""
    seq_ftr = create_2_part_seqfeature(152829, 153401, 1, 0, 4, 1, "CDS")
    return seq_ftr




# Alice CDS 124
# Example 2-part compound feature on top strand.

alice_cds_124_translation = (
    "MKDEMTTSDVPADPAIDPDLAPPEPRRVVGELVETEPQEHEDPEVTELTDEERSSFVSLLT"
    "CGKHSKKITVMGHPVVIQTLKTGDEMRVGLFTKKYLESQMGFQRAYQVAVCAAGIREIQGK"
    "PLFRELREVTDEDEIFDKNVEAVMELYPIVITQIYQAIMDLEREYAQLAVKLGKTVRLDAS"
    "TELEIRLAYKQGLLTQPSLNRYQRWALRYAIFMDRRLQLQDTEDMLQRQTWYLEPKRYHDL"
    "FLAGAFEPEPIAVAGRDMEEVVDDLDEVDAYFARLEGSQSMSGAQLFAALDEPDEEGWM"
    )

def get_alice_cds_124_draft_data_in_db():
    """Returns mock data from gene table for 'draft' Alice CDS 124."""
    dict = {
        "GeneID": "Alice_CDS_1",
        "PhageID": "Alice",
        "Start": 70374,
        "Stop": 71285,
        "Parts": 2,
        "Length": (len(alice_cds_124_translation) * 3) + 3,
        "Name": "124",
        "Translation": alice_cds_124_translation,
        "Orientation": "F",
        "Notes": "",
        "LocusTag": None,
        "DomainStatus": 0,
        "PhamID": "NULL"
        }
    return dict

def get_alice_cds_124_qualifier_dict():
    """Returns data to construct Alice CDS 124 flat file feature."""
    dict = OrderedDict(
        [("gene", ["124"]),
         ("locus_tag", ["ALICE_124"]),
         ("ribosomal_slippage", [""]),
         ("note", ["gp124"]),
         ("codon_start", ["1"]),
         ("transl_table", ["11"]),
         ("product", [""]),
         ("protein_id", ["AEJ94379.1"]),
         ("translation", [])])
    return dict

def get_alice_cds_124_seqfeature():
    """Constructs Alice CDS 124 flat file feature."""
    seq_ftr = create_2_part_seqfeature(70374, 70902, 1, 70901, 71285, 1, "CDS")
    return seq_ftr




# Alice tmrna_169
# Example tmRNA feature.
def get_alice_tmrna_169_draft_data_in_db():
    """Returns mock data from trna table for 'draft' Alice tmRNA 169."""
    dict = {
        "GeneID": "Alice_TMRNA_1",
        "PhageID": "Alice",
        "Start": 95923,
        "Stop": 96358,
        "Length": 435,
        "Name": "169",
        "Orientation": "F",
        "Note": "Tag peptide: ATDTDATVTDAEIEAFFAEEAAALV*",
        "LocusTag": "NULL",
        "PeptideTag": "ATDTDATVTDAEIEAFFAEEAAALV*"
        }
    return dict


def get_alice_tmrna_169_qualifier_dict():
    """Returns data to construct Alice tmRNA 169 flat file feature."""
    dict = OrderedDict(
            [('gene', ['169']),
             ('locus_tag', ['ALICE_169']),
             ('note', ['Tag peptide: ATDTDATVTDAEIEAFFAEEAAALV*'])])
    return dict

def get_alice_tmrna_169():
    """Constructs Alice tmRNA 169 flat file feature."""
    seq_ftr = create_1_part_seqfeature(95923, 96358, 1, "tmRNA")
    return seq_ftr




# Alice trna_170
# Example tRNA feature.
def get_alice_trna_170_draft_data_in_db():
    """Returns mock data from trna table for 'draft' Alice tRNA 170."""
    dict = {
        "GeneID": "Alice_TRNA_1",
        "PhageID": "Alice",
        "Start": 96431,
        "Stop": 96507,
        "Length": 76,
        "Name": "170",
        "Orientation": "F",
        "Note": "tRNA-Gln(ttg)",
        "LocusTag": "NULL",
        "AminoAcid": "Gln",
        "Anticodon": "ttg",
        "Structure": "NULL",
        "Source": "NULL",
        }
    return dict


def get_alice_trna_170_qualifier_dict():
    """Returns data to construct Alice tRNA 170 flat file feature."""
    dict = OrderedDict(
            [('gene', ['170']),
             ('locus_tag', ['ALICE_170']),
             ('product', ['tRNA-Gln']),
             ('note', ['tRNA-Gln(ttg)'])])
    return dict

def get_alice_trna_170():
    """Constructs Alice tRNA 170 flat file feature."""
    seq_ftr = create_1_part_seqfeature(96431, 96507, 1, "tRNA")
    return seq_ftr




# Alice source_1
# Example source feature.

def get_alice_source_1_qualifiers():
    """Returns data to construct Alice source flat file feature."""
    dict = OrderedDict(
            [("organism", ["Mycobacterium phage Alice_Draft"]),
             ("mol_type", ["genomic DNA"]),
             ("isolation_source", ["soil"]),
             ("db_xref", ["taxon:1034128"]),
             ("lab_host", ["Mycobacterium smegmatis mc2 155"]),
             ("country", ["USA: Aledo, TX"]),
             ("lat_lon", ["31.69 N 97.65 W"]),
             ("collection_date", ["01-Oct-2009"]),
             ("collected_by", ["C. Manley"]),
             ("identified_by", ["C. Manley"])])
    return dict

def get_alice_source_1():
    """Constructs Alice source flat file feature."""
    seq_ftr = create_1_part_seqfeature(0, 153401, 1, "source")
    return seq_ftr




# Authors

author_string_1 = "Hatfull,G.F."
author_string_2 = (
    "Alferez,G.I., Bryan,W.J., Byington,E.L., Contreras,T.D., "
    "Evans,C.R., Griffin,J.A., Jalal,M.D., Lindsey,C.B., "
    "Manley,C.M., Mitchell,A.M., O'Hara,J.M., Onoh,U.M., "
    "Padilla,E., Penrod,L.C., Regalado,M.S., Reis,K.E., "
    "Ruprecht,A.M., Slater,A.E., Staton,A.C., Tovar,I.G., "
    "Turek,A.J., Utech,N.T., Simon,S.E., Hughes,L.E., "
    "Benjamin,R.C., Serrano,M.G., Lee,V., Hendricks,S.L., "
    "Sheth,N.U., Buck,G.A., Bradley,K.W., Khaja,R., "
    "Lewis,M.F., Barker,L.P., Jordan,T.C., Russell,D.A., "
    "Leuba,K.D., Fritz,M.J., Bowman,C.A., Pope,W.H., "
    "Jacobs-Sera,D., Hendrix,R.W. and Hatfull,G.F.")




###Trixie genome
# Data is not accurate, but primarily serves as simple test data.

def get_trixie_phage_data():
    """Mock phage table data for Trixie.

    Data is not accurate, but serves as simple test data."""
    dict = {
        "PhageID": "Trixie",
        "Accession": "BCD456",
        "Name": "Trixie",
        "HostGenus": "Gordonia",
        "Sequence": "GGGGGGGGGGGGGGGGGGGG",
        "Length": 20,
        "GC": 1,
        "Status": "final",
        "DateLastModified": datetime.strptime('1/1/2000', '%m/%d/%Y'),
        "RetrieveRecord": 1,
        "AnnotationAuthor": 1,
        "Cluster": "A",
        "Subcluster": "A3",
        "Notes": "NULL"
        }
    return dict

def get_trixie_gene_data():
    """Mock gene table data for Trixie.

    Data is not accurate, but serves as simple test data."""
    dict = {
        "GeneID": "TRIXIE_0001",
        "PhageID": "Trixie",
        "Start": 100,
        "Stop": 1100,
        "Parts": 1,
        "Length": 1000,
        "Name": "1",
        "Translation": "ACTGC",
        "Orientation": "F",
        "Notes": "int",
        "LocusTag": "SEA_TRIXIE_0001",
        "DomainStatus": 0,
        "PhamID": "NULL"
        }
    return dict

def get_trixie_gene_domain_data():
    """Mock gene_domain table data for Trixie.

    Data is not accurate, but serves as simple test data."""
    dict = {
        "GeneID": "TRIXIE_0001",
        "HitID": "gnl|CDD|334841",
        "Expect": 1.78531e-11,
        "QueryStart": 33,
        "QueryEnd": 115
        }
    return dict

def get_trixie_domain_data():
    """Mock domain table data for Trixie.

    Data is not accurate, but serves as simple test data."""
    dict = {
        "HitID": "gnl|CDD|334841",
        "DomainID": "pfam02195",
        "Name": "ParBc",
        "Description": "ParB-like nuclease domain"
        }
    return dict

def get_trixie_trna_data():
    """Mock trna table data for Trixie.

    Data is not accurate, but serves as simple test data."""
    dict = {
        "GeneID": "TRIXIE_0001",
        "PhageID": "Trixie",
        "Start": 100,
        "Stop": 1100,
        "Length": 1000,
        "Name": "1",
        "Orientation": "F",
        "Note": "misc",
        "LocusTag": "SEA_TRIXIE_0001",
        "AminoAcid": "Ala",
        "Anticodon": "AAA",
        "Structure": "AAAAAAAA",
        "Source": "aragorn",
        }
    return dict

def get_trixie_tmrna_data():
    """Mock tmRNA table data for Trixie.

    Data is not accurate, but serves as simple test data."""
    dict = {
        "GeneID": "TRIXIE_0001",
        "PhageID": "Trixie",
        "Start": 100,
        "Stop": 1100,
        "Length": 1000,
        "Name": "1",
        "Orientation": "F",
        "Note": "misc",
        "LocusTag": "SEA_TRIXIE_0001",
        "PeptideTag": "random"
        }
    return dict


###Lifes genome
# Lifes_Draft CDS 122
# Example 2-part compound feature on bottom strand that wraps around genome termini.

def get_lifes_cds_122_seqfeature():
    """Constructs Lifes_Draft CDS 122 flat file feature."""
    # Below is the output when BioPython parses this feature.
    # SeqFeature(
    #     CompoundLocation(
    #         [FeatureLocation(
    #             ExactPosition(0),
    #             ExactPosition(9),
    #             strand=-1
    #             ),
    #          FeatureLocation(
    #             ExactPosition(58743),
    #             ExactPosition(59253),
    #             strand=-1
    #             )
    #          ], 'join'
    #      ),
    #      type='CDS',
    #      location_operator='join')
    seq_ftr = create_2_part_seqfeature(0, 9, -1, 58743, 59253, -1, "CDS")
    return seq_ftr
