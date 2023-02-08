import os

# create log directory - needed for cluster slurm execution
os.makedirs("log", exist_ok=True)


configfile: "config/config.yml"


include: "rules/downloads.smk"
include: "rules/pangraph.smk"
include: "rules/distances.smk"


localrules:
    download_gbk,
    PG_all,
