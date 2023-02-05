# Workflow description

## download.smk

Contains rules for downloading data given an accession number `{acc}`.

```mermaid
flowchart TD
    A("{acc}") --> |"download_gbk"| B("data/gbk/{acc}.gbk") --> |"gbk_to_fa"| C("data/fa/{acc}.fa")
```

## pangraph.smk

Rules to build, polish and export a pangenome graph given a `{dset}` (collection of accession numbers) and a kernel option `{opt}`.

```mermaid
flowchart TD
    A("data/fa/{acc}.fa from {dset}") --> |PG_build| B("pangraph/{dset}/{opt}.json") 
    B --> |PG_polish| C("pangraph/{dset}/{opt}-polished.json")
    C --> |PG_export| D("pangraph/{dset}/export/{opt}")
    C --> |PG_summary_fig| E["figs/{dset}/pangraph/{opt}_summary.pdf"]
```

**Description**:
- `pangraph/{dset}/export/{opt}` : folder containing `.gfa` export of the polished pangraph.
- `figs/{dset}/pangraph/{opt}_summary.pdf` : summary figure with distribution of block frequency/length.
