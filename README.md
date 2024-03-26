# Pangenome evolution

This repository contains scripts for the analysis of microbial pangenome evolution.

## setup

Execution requires a valid installation of conda, mamba and snakemake (v7.32.4).
For pangenome graph creation, the `pangraph` command must be available in path.

### ncbi api key

To facilitate download of genbank records from ncbi, your personal api key can be saved in `config/ncbi_api_key.txt`. It will be automatically used when downloading the data.

### for snakemake v8+

For `slurm` execution, install the [slurm executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html):
```bash
pip install snakemake-executor-plugin-slurm
```
and execute with:
```bash
snakemake --profile slurm <rule>
```

