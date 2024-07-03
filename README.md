# Pangenome evolution

This repository contains a pipeline for the analysis of pangenome evolution

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
