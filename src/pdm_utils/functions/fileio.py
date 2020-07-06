from pathlib import Path

def write_fasta(ids_seqs, infile_path):
    """
    Writes the input genes to the indicated file in FASTA multiple
    sequence format (unaligned).
    :param id_seqs: the ids and sequences to be written to file 
    :type genes: dict
    :param infile_path: the path of the file to write the genes to
    :type infile: Path
    :type infile: str
    """
    if isinstance(infile, str):
        file_handle = open(infile_path, "w")
    elif isinstance(infile_path, Path):
        file_handle = infile.open(mode="w")
    else:
        raise TypeError("File path type not supported.")

    for id, seq in ids_seqs.items():
        split_seq = textwrap.wrap(seq, 60)
        wrapped_seq = "\n".join(split_seq)
        file_handle.write(f">{id}\n{wrapped_seq}\n")

    file_handle.close()
