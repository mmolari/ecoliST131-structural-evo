from Bio import SeqIO
import pandas as pd
import json


def gbk_to_loci_df(gbk_file):
    gb = SeqIO.read(gbk_file, "genbank")

    iso = gb.annotations["source"]

    len_genome = len(gb.seq)

    # create a dataframe of loci
    loci = []
    for f in gb.features:
        if f.type == "CDS":
            loc = f.location
            # check if compount positions
            parts = loc.parts
            if len(parts) > 1:

                # check that start < end for every part
                for p in parts:
                    assert (
                        p.start < p.end
                    ), f"start: {p.start} | end: {p.end} | {f.qualifiers['locus_tag'][0]} | {iso} | {gbk_file}"

                # pick minimum and maximum for start/end
                start = min([p.start for p in parts])
                end = max([p.end for p in parts])

                # if these span more than 1 Mbp then there is something wrong (wrap around the genome)
                if end - start > 1e6:
                    print(
                        f"wraparound for {gbk_file} | {f.qualifiers['locus_tag'][0]} | {parts}"
                    )
                    start = min([p.start for p in parts if p.start > len_genome / 2])
                    end = max([p.end for p in parts if p.end < len_genome / 2])
                    print(f"start: {start} | end: {end}")
            else:
                start = loc.start
                end = loc.end
            loci.append(
                {
                    "locus": f.qualifiers["locus_tag"][0],
                    "start": start,
                    "end": end,
                    "strand": loc.strand,
                    "length": len(f),
                }
            )
    loci = pd.DataFrame(loci)
    return loci


def load_genecluster_json(fname):
    with open(fname, "r") as f:
        jc = pd.DataFrame(json.load(f))
    jc["geneId"] = jc["geneId"].astype(int)
    jc["divers"] = jc["divers"].astype(float)
    jc["count"] = jc["count"].astype(int)
    jc["geneLen"] = jc["geneLen"].astype(int)
    jc["event"] = jc["event"].astype(int)
    jc.set_index("geneId", inplace=True)
    return jc
