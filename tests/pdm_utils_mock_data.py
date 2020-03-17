"""Mock Alice genome data for building mock databases and mock flat files."""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
# TODO utilize BeforePosition and Reference
from Bio.SeqFeature import ExactPosition, BeforePosition, Reference
from collections import OrderedDict
from datetime import datetime
from pathlib import Path


# Alice ("test_flat_file_10.gb"),
unittest_file = Path(__file__)
test_dir = unittest_file.parent
test_file_dir = Path(test_dir, "test_files")
base_flat_file = Path("test_flat_file_10.gb")
base_flat_file_path = Path(test_file_dir, base_flat_file)






# Format of the date the script imports into the database.
current_date = datetime.today().replace(hour=0, minute=0,
                                        second=0, microsecond=0)


def get_seq(filepath):
    """Get genome sequence from a flat file so that it can be added
    to a MySQL database for testing."""
    seqrecord = SeqIO.read(filepath, "genbank")
    return seqrecord.seq



alice_cds_252_translation = (
    "MFDNHPLMSTLITSVRDLNRLANGVIIRNRCEKCDTSAGPDQGV"
    "HFVKTPLGWLYSDPKGKPTTWELFPSDAIHLPVQVIHRVSESEDVSEITENSPHTRKD"
    "DSEATEGAAPSFEPLYPSSHKTERFTVQDTPDGLGAYVFDFGGDAFGAQTCADALAKV"
    "TGKTWYVMHKTVVENTGIYSSMVRHEEESSS")


def get_alice_cds_252_draft_data_in_db(translation=alice_cds_252_translation):
    """Returns mock data from gene table for 'draft' Alice CDS 252."""
    dict = {
        "GeneID": "Alice_CDS_4",
        "PhageID": "Alice",
        "Start": 152829,
        "Stop": 4,
        "Parts": 2,
        "Length": 191,
        "Name": "252",
        "Translation": translation,
        "Orientation": "F",
        "Notes": "",
        "LocusTag": None
        }
    return dict


# TODO this may no longer be needed.
def get_alice_cds_252_final_data_in_db():
    """Returns mock data from gene table for 'final' Alice CDS 252."""
    dict = get_alice_cds_252_draft_data_in_db()
    dict["LocusTag"] = "ALICE_252"
    return dict


alice_cds_124_translation = (
    "MKDEMTTSDVPADPAIDPDLAPPEPRRVVGELVETEPQEHEDPEVTELTDEERSSFVSLLT"
    "CGKHSKKITVMGHPVVIQTLKTGDEMRVGLFTKKYLESQMGFQRAYQVAVCAAGIREIQGK"
    "PLFRELREVTDEDEIFDKNVEAVMELYPIVITQIYQAIMDLEREYAQLAVKLGKTVRLDAS"
    "TELEIRLAYKQGLLTQPSLNRYQRWALRYAIFMDRRLQLQDTEDMLQRQTWYLEPKRYHDL"
    "FLAGAFEPEPIAVAGRDMEEVVDDLDEVDAYFARLEGSQSMSGAQLFAALDEPDEEGWM"
    )


def get_alice_cds_124_draft_data_in_db(translation=alice_cds_124_translation):
    """Returns mock data from gene table for 'draft' Alice CDS 124."""
    dict = {
        "GeneID": "Alice_CDS_1",
        "PhageID": "Alice",
        "Start": 70374,
        "Stop": 71285,
        "Parts": 2,
        "Length": 303,
        "Name": "124",
        "Translation": translation,
        "Orientation": "F",
        "Notes": "",
        "LocusTag": None
        }
    return dict


# TODO this may no longer be needed.
def get_alice_cds_124_final_data_in_db():
    """Returns mock data from gene table for 'final' Alice CDS 124."""
    dict = get_alice_cds_124_draft_data_in_db()
    dict["LocusTag"] = "ALICE_124"
    dict["Notes"] = "tail assembly chaperone"
    return dict


alice_cds_139_translation = (
    "MKKLVTGGLVAAGITLGLISAPSASAEYMTICPSQVSAVVTANTSCGFADNV"
    "FRGFYRQSGWDPLAYSPATGKVYRMHCAPATTTNWGEAKRCWGVGYGGDLLVVYID"
    )


def get_alice_cds_139_draft_data_in_db(translation=alice_cds_139_translation):
    """Returns mock data from gene table for 'draft' Alice CDS 139."""
    dict = {
        "GeneID": "Alice_CDS_2",
        "PhageID": "Alice",
        "Start": 88120,
        "Stop": 88447,
        "Parts": 1,
        "Length": 108,
        "Name": "139",
        "Translation": translation,
        "Orientation": "R",
        "Notes": "",
        "LocusTag": None
        }
    return dict

# TODO this may no longer be needed.
def get_alice_cds_139_final_data_in_db():
    """Returns mock data from gene table for 'final' Alice CDS 139."""
    dict = get_alice_cds_139_draft_data_in_db()
    dict["LocusTag"] = "ALICE_139"
    return dict


alice_cds_193_translation = (
    "MGNNRPTTRTLPTGEKATHHPDGRLVLKPKNSLADALSGQLDAQQEASQALTEALV"
    "QTAVTKREAQADGNPVPEDKRVF"
    )


def get_alice_cds_193_draft_data_in_db(translation=alice_cds_193_translation):
    """Returns mock data from gene table for 'draft' Alice CDS 193."""
    dict = {
        "GeneID": "Alice_CDS_3",
        "PhageID": "Alice",
        "Start": 110297,
        "Stop": 110537,
        "Parts": 1,
        "Length": 79,
        "Name": "193",
        "Translation": translation,
        "Orientation": "F",
        "Notes": "",
        "LocusTag": None
        }
    return dict

# TODO this may no longer be needed
def get_alice_cds_193_final_data_in_db():
    """Returns mock data from gene table for 'final' Alice CDS 193."""
    dict = get_alice_cds_193_draft_data_in_db()
    dict["LocusTag"] = "ALICE_193"
    return dict


# TODO not sure if this belongs in this module.
# Tuples of the four Alice CDS features used.
# alice_cds_252_coords = (152829, 4)
# alice_cds_124_coords = (70374, 71285)
# alice_cds_139_coords = (88120, 88447)
# alice_cds_193_coords = (110297, 110537)


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
    seq_ftr = SeqFeature(
                CompoundLocation(
                    [FeatureLocation(
                        ExactPosition(152829),
                        ExactPosition(153401),
                        strand=1),
                    FeatureLocation(
                        ExactPosition(0),
                        ExactPosition(4),
                        strand=1)],
                    "join"),
                type="CDS",
                location_operator="join")
    return seq_ftr


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
    seq_ftr = SeqFeature(
                CompoundLocation(
                    [FeatureLocation(
                        ExactPosition(70374),
                        ExactPosition(70902),
                        strand=1),
                    FeatureLocation(
                        ExactPosition(70901),
                        ExactPosition(71285),
                        strand=1)],
                    "join"),
                type="CDS",
                location_operator="join")
    return seq_ftr


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
    seq_ftr = SeqFeature(
                FeatureLocation(
                    ExactPosition(88120),
                    ExactPosition(88447),
                    strand=-1),
                type="CDS")
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


def get_alice_cds_193_seqfeature():
    """Constructs Alice CDS 193 flat file feature."""
    seq_ftr = SeqFeature(
                FeatureLocation(
                    ExactPosition(110297),
                    ExactPosition(110537),
                    strand=1),
                type="CDS")
    return seq_ftr


def get_alice_tmrna_169_qualifier_dict():
    """Returns data to construct Alice tmRNA 169 flat file feature."""
    dict = OrderedDict(
            [('gene', ['169']),
             ('locus_tag', ['ALICE_169'])])
    return dict


def get_alice_tmrna_169():
    """Constructs Alice tmRNA 169 flat file feature."""
    seq_ftr = SeqFeature(
                FeatureLocation(
                    ExactPosition(95923),
                    ExactPosition(96358),
                    strand=1),
                type="tmRNA")
    return seq_ftr


def get_alice_trna_170_qualifier_dict():
    """Returns data to construct Alice tRNA 170 flat file feature."""
    dict = OrderedDict(
            [('gene', ['170']),
             ('locus_tag', ['ALICE_170']),
             ('product', ['tRNA-Gln']),
             ('note', ['tRNA-Gln (ttg)'])])
    return dict


def get_alice_trna_170():
    """Constructs Alice tRNA 170 flat file feature."""
    seq_ftr = SeqFeature(
                FeatureLocation(
                    ExactPosition(96431),
                    ExactPosition(96507),
                    strand=1),
                type="tRNA")
    return seq_ftr


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
    seq_ftr = SeqFeature(
                FeatureLocation(
                    ExactPosition(0),
                    ExactPosition(153401),
                    strand=1),
                type="source")
    return seq_ftr


# Authorship info.
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




def get_alice_genome_draft_data_in_db(seq=True):
    """Returns mock data from phage table for 'draft' Alice genome."""
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















#
