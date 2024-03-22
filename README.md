# Pangenome evolution

This repository contains scripts for the analysis of microbial pangenome evolution.

## setup

Execution requires a valid installation of conda, mamba and snakemake.
For pangenome graph creation, the `pangraph` command must be available in path.

### for snakemake v8+

For `slurm` execution, install the [slurm executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html):
```bash
pip install snakemake-executor-plugin-slurm
```
and execute with:
```bash
snakemake --profile slurm <rule>
```
