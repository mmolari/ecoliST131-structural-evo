# Structural genome evolution in _E. coli_ ST131

This repository contains a snakemake pipeline for the analysis of structural genomic evolution of _E.coli_ ST131 presented in [our paper](#citation).

The **dataset** consists of complete _E. coli_ ST131 genomes available on [RefSeq](https://www.ncbi.nlm.nih.gov/datasets/genome/). Accession numbers and metadata for the considered strains can be found in [the datasets folder](config/datasets/ST131_ABC).

In short, the **pipeline** uses [pangraph](https://github.com/neherlab/pangraph) to build a pangenome graph representation for the chromosomes of all of the considered strains. It then extracts all regions of structural variations, assigns MGEs and defense systems to each of these regions, and detect events that can be parsimoniously interpreted as simple gain or loss of sequence. See [this note for an overview of the pipeline](notes/workflow.md).

The pipeline produces as **output** a `results` folder, containing processed data such as the pangenome graph and the junction graphs, and a `figs` folder, containing amongst other the main figures of [the paper](#citation).

## setup

- Execution requires a valid installation of [conda](https://conda.io/projects/conda), [mamba](https://mamba.readthedocs.io) and [snakemake](https://snakemake.readthedocs.io) (v7.32.4).
- For pangenome graph creation, the [pangraph](https://github.com/neherlab/pangraph) command must be available in path, see [pangraph documentation](https://neherlab.github.io/pangraph/#Installation) for installation instructions.
- optionally, to facilitate download of genbank records from ncbi, your personal api key can be saved in `config/ncbi_api_key.txt`. It will be automatically used when downloading the data.

## execution

to execute the pipeline locally, it is sufficient to run:
```sh
snakemake --use-conda --cores 1 all
```
You can replace `1` with the desired number of cores.

Give the high number of jobs and the memory and time requirements we advise executing on cluster. Execution using the SLURM workload manager is already set up and the pipeline can be executed with:
```sh
snakemake --profile cluster all
```

## citation

**Evolutionary dynamics of genome structure and content among closely related bacteria**

Marco Molari, Liam P. Shaw and Richard A. Neher, _biorxiv_ (2024)