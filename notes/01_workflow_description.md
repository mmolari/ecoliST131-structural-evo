# Workflow description

## download.smk

Contains rules for downloading data given an accession number `{acc}`.

```mermaid
flowchart TD
    A("{acc}") --> |"download_gbk"| B("data/gbk/{acc}.gbk") --> |"gbk_to_fa"| C("data/fa/{acc}.fa")
```

## pangraph.smk

Rules to build, polish and export a pangenome graph given a `{dset}` (collection of accession numbers) and a kernel option `{opt}`. Results are saved either in the `figs/{dset}` or `results/{dset}` subfolders

```mermaid
flowchart TD
    A("data/fa/{acc}.fa from {dset}") --> |PG_build| B("{dset}/pangraph/{opt}.json") 
    B --> |PG_polish| C("pangraph/{opt}-polished.json")
    C --> |PG_export| D("pangraph/export/{opt}")
    C --> |PG_block_distr_fig| E["pangraph/{opt}_block_distr.pdf"]
    C --> |PG_reduced_corealignment| F("pangraph/{opt}-alignment/corealignment{.fa,_info.json}")
```

**Description**:
- `pangraph/export/{opt}` : folder containing `.gfa` export of the polished pangraph.
- `pangraph/{opt}_block_distr.pdf` : figure with distribution of block frequency/length.
- `pangraph/{opt}-alignment/corealignment` : reduced core alignment, and info file with number of sites having gaps / being consensus.