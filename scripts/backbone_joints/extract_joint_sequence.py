from Bio import SeqIO
import json
import pathlib

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        given the dictionary of edge positions and the genomes fasta files,
        it extracts from each genome the sequence corresponding to the desired
        interval and saves it in a fasta file.
        """
    )
    parser.add_argument("--edge_pos", type=str)
    parser.add_argument("--edge", type=str)
    parser.add_argument("--out_fa", type=str)
    parser.add_argument("--genomes", type=str, nargs="+")
    return parser.parse_args()


# create dictionary isolate -> fasta file for input genomes:
def genomes_file_dict(fnames):
    file_dict = {}
    for fname in fnames:
        iso = pathlib.Path(fname).stem
        file_dict[iso] = fname
    return file_dict


if __name__ == "__main__":
    args = parse_args()

    # input genome files
    genomes_fdict = genomes_file_dict(args.genomes)

    # load edge positions
    with open(args.edge_pos, "r") as f:
        all_edge_pos = json.load(f)
    edge_pos = all_edge_pos[args.edge]

    records = []
    for iso, pos in edge_pos.items():
        beg, _, _, end, strand = pos

        # load genome
        with open(genomes_fdict[iso], "r") as f:
            seq = SeqIO.read(f, "fasta")

        # if it wraps around the genome
        if beg < end:
            sseq = seq[beg - 1 : end]
        else:
            sseq = seq[beg - 1 :] + seq[:end]

        # reverse-complement if necessary
        if not strand:
            sseq = sseq.reverse_complement()

        # append record
        sseq.id = f"{iso}"
        sseq.description = ""
        records.append(sseq)

    with open(args.out_fa, "w") as f:
        SeqIO.write(records, f, "fasta")
