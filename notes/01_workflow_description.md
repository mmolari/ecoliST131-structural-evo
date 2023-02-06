# Workflow description

## download.smk

Contains rules for downloading data given an accession number `{acc}`.

```mermaid
flowchart TD
    A("{acc}") --> |"download_gbk"| B("data/gbk/{acc}.gbk") --> |"gbk_to_fa"| C("data/fa/{acc}.fa")
```

## pangraph.smk

Rules to build, polish and export a pangenome graph given a `{dset}` (collection of accession numbers) and a kernel option `{opt}`. Results are saved either in the `figs/{dset}/pangraph` or `results/{dset}/pangraph` subfolders, depending on their type.

```mermaid
flowchart TD
    A("data/fa/{acc}.fa from {dset}") --> |PG_build| B("{opt}.json") 
    B --> |PG_polish| C("{opt}-polished.json")
    C --> |PG_export| D("export/{opt}")
    C --> |PG_block_distr_fig| E["{opt}_block_distr.pdf"]
    C --> |PG_reduced_corealignment| F("{opt}-alignment/corealignment{.fa,_info.json}")
    F --> |PG_coregenome_tree| G("{opt}-coretree.nwk")
```

**Description**:
- `export/{opt}` : folder containing `.gfa` export of the polished pangraph.
- `{opt}_block_distr.pdf` : figure with distribution of block frequency/length.
- `{opt}-alignment/corealignment{.fa,_info.json}` : reduced core alignment, and info file with number of sites having gaps / being consensus.
- `{opt}-coretree.nwk` : core genome tree build from the reduced alignment (rescaled with information on the number of consensus vs mutated sites).

## distances.smk

Rules to estimate evolutionary distances between strains. Results are saved in `results/{dset}/distances`.

```mermaid
flowchart TD
    A("{opt}-alignment/corealignment{.fa,_info.json}") --> |DST_corealignment| B("{opt}-coredivergence.csv")
    C("data/fa/{acc}.fa from {dset}") --> |DST_mash| D("mash_dist.csv")
```

**Description**
- `{opt}-coredivergence.csv` : core genome divergence, evaluated from the restricted core genome alignment, and rescaled with alignment info to the full core genome.
- `mash_dist.csv` : mash distance.