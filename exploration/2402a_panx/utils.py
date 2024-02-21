from Bio import SeqIO
import pandas as pd

def gbk_to_loci_df(gbk_file):
    gb = SeqIO.read(gbk_file, "genbank")

    iso = gb.annotations["source"]

    # create a dataframe of loci
    loci = []
    for f in gb.features:
        if f.type == "CDS":
            loc = f.location
            # check if compount positions
            if len(loc.parts) > 1:
                # print(f.qualifiers["locus_tag"][0], "compound")
                # print(loc.parts)
                start = loc.parts[0].start
                end = loc.parts[-1].end
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