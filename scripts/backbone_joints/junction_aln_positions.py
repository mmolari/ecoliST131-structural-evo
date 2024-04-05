import numpy as np
import pandas as pd
import pypangraph as pp
import utils as ut

import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--guide_strain")
    parser.add_argument("--thr_len", type=int)
    parser.add_argument("--genome_lengths")
    parser.add_argument("--pangraph")
    parser.add_argument("--out_SNPs")
    parser.add_argument("--out_junct_pos")
    parser.add_argument("--out_aln_pos")
    return parser.parse_args()


def load_guide_genome_length(args):
    genome_lengths = pd.read_csv(args.genome_lengths, index_col=0)
    guide_L = genome_lengths.loc[args.guide_strain]["length"]
    return guide_L


def load_guide_strain_path(pan, bdf, guide_strain):
    p = pan.paths[guide_strain]
    blocks = p.block_ids
    strands = p.block_strands
    positions = p.block_positions[:-1]

    # rotate to longest core-block
    lcb = bdf[bdf["core"]]["len"].idxmax()
    lcb_idx = np.argwhere(blocks == lcb)[0][0]
    blocks = np.roll(blocks, -lcb_idx)
    strands = np.roll(strands, -lcb_idx)
    positions = np.roll(positions, -lcb_idx)

    return blocks, strands, positions


def ungapped_aln(aln_dict, strains):

    L = len(aln_dict[strains[0]])
    N = len(strains)
    M = np.empty((N, L), dtype=str)

    for n, s in enumerate(strains):
        M[n] = list(aln_dict[s])

    # remove gaps
    has_gap = np.any(M == "-", axis=0)
    M = M[:, ~has_gap]

    # remove extra characters
    is_extra = ~np.all(np.isin(M, list("acgtACGT")), axis=0)
    M = M[:, ~is_extra]

    return M


if __name__ == "__main__":

    args = parse_args()

    # load data
    guide_L = load_guide_genome_length(args)
    pan = pp.Pangraph.load_json(args.pangraph)
    strains = pan.strains()
    bdf = pan.to_blockstats_df()

    # load guide path
    blocks, strands, positions = load_guide_strain_path(pan, bdf, args.guide_strain)

    snps_info = {}
    junct_pos = {}
    aln_pos = {}
    core_left = None
    core_right = None
    L_current = 0
    for nbl, (b, s, pos) in enumerate(zip(blocks, strands, positions)):
        b_len = bdf.loc[b]["len"]
        if (not bdf.loc[b]["core"]) or (b_len < args.thr_len):
            continue
        core_strand = pan.paths[args.guide_strain]

        fr = nbl / len(blocks)
        print(f"completed {fr*100:.2f} % |{int(fr*10)*'#':<10}|", end="\r")
        # left and right positions on guide strain

        cbl = pan.blocks[b]
        alns, occs = cbl.alignment.generate_alignments()

        # alignment dictionary
        aln_dict = {occ[0]: aln for aln, occ in zip(alns, occs)}

        # ungapped alignment
        M = ungapped_aln(aln_dict, strains)

        # reduced length
        L_aln = M.shape[1]

        # find non-consensus columns
        is_consensus = np.all(M == M[0], axis=0)
        SNPs_pos = np.argwhere(~is_consensus).flatten() + L_current
        M_red = M[:, ~is_consensus]

        # save SNPs
        for p, M_red_col in zip(SNPs_pos, M_red.T):
            snps_info[p] = M_red_col

        # assign junction breakpoint
        if core_left is None:
            core_left = ut.Node(b, s)
        else:
            core_right = ut.Node(b, s)
            j = ut.Edge(core_left, core_right)
            j = j.to_str_id()
            junct_pos[j] = L_current + L_aln
            core_left = core_right

        # assign alignment positions
        aln_pos[b] = {
            "start": L_current,
            "end": L_current + L_aln,
            "len": L_aln,
            "guide_start": pos,
            "guide_end": positions[(nbl + 1) % len(blocks)],
        }
        aln_pos[b]["guide_len"] = (
            aln_pos[b]["guide_end"] - aln_pos[b]["guide_start"]
        ) % guide_L

        L_current += L_aln

    # assign the last one
    core_right = ut.Node(blocks[0], strands[0])
    j = ut.Edge(core_left, core_right)
    j = j.to_str_id()
    junct_pos[j] = L_current

    # convert to dataframes
    junct_pos = pd.DataFrame(junct_pos.items(), columns=["junction", "position"])
    junct_pos = junct_pos.set_index("junction")
    snps_info = pd.DataFrame(snps_info).T
    snps_info.columns = strains
    aln_pos = pd.DataFrame(aln_pos).T

    # save output
    junct_pos.to_csv(args.out_junct_pos)
    snps_info.to_csv(args.out_SNPs)
    aln_pos.to_csv(args.out_aln_pos)
