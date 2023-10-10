import os

# create log directory - needed for cluster slurm execution
os.makedirs("log", exist_ok=True)


configfile: "config/config.yml"


# load config file entries
dsets_config = config["datasets"]
dset_names = list(dsets_config.keys())

# read accession numbers from dataset files
dset_chrom_accnums = {}
for dset_name, dset_info in dsets_config.items():
    fname = dset_info["chromosomes"]
    with open(fname, "r") as f:
        acc_nums = f.readlines()
    acc_nums = [an.strip() for an in acc_nums]
    acc_nums = [an for an in acc_nums if len(an) > 0]
    dset_chrom_accnums[dset_name] = acc_nums

# load accession numbers of excluded isolates
excluded = {k: [] for k in dset_names}
for dset_name, dset_info in dsets_config.items():
    fname = dset_info["excluded"]
    with open(fname, "r") as f:
        acc_nums = f.readlines()
    acc_nums = [an.strip() for an in acc_nums]
    acc_nums = [an for an in acc_nums if len(an) > 0]
    excluded[dset_name] = acc_nums

# pangraph kernel options
kernel_opts = list(config["pangraph"]["kernel-options"].keys())


wildcard_constraints:
    opt=f"({'|'.join(kernel_opts)})",
    dset=f"({'|'.join(dset_names)})",


include: "rules/downloads.smk"
include: "rules/pangraph.smk"
include: "rules/resistance.smk"
include: "rules/distances.smk"
include: "rules/assembly_qc.smk"
include: "rules/backbone_joints.smk"
include: "rules/figs.smk"


localrules:
    download_gbk,
    QC_busco_download,
