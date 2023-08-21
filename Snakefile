import os

# create log directory - needed for cluster slurm execution
os.makedirs("log", exist_ok=True)


configfile: "config/config.yml"


include: "rules/downloads.smk"
include: "rules/pangraph.smk"
include: "rules/distances.smk"
include: "rules/assembly_qc.smk"
include: "rules/backbone_joints.smk"


localrules:
    download_gbk,
    QC_busco_download,
