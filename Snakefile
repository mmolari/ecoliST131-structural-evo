import os
import json

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

# load plasmids
# {dset_name: {chromosome_accnum: [plasmid_accnum_1, plasmid_accnum_2 ...]}}
plasmid_dsets = {}
for dset_name in dsets_config:
    if "plasmids" in dsets_config[dset_name]:
        json_file = dsets_config[dset_name]["plasmids"]
        with open(json_file, "r") as f:
            plasmid_dsets[dset_name] = json.load(f)

# load plasmid accnums:
# {dset_name: [accnum1, accnum2, ...]}
plasmid_accnums = {}
for dset_name, pl_dset in plasmid_dsets.items():
    plasmid_accnums[dset_name] = sum(pl_dset.values(), [])

# pangraph kernel options
kernel_opts = list(config["pangraph"]["kernel-options"].keys())


wildcard_constraints:
    opt="(" + "|".join(kernel_opts) + ")",
    dset="(" + "|".join(dset_names) + ")",
    acc=r"[^/]+",


include: "rules/downloads.smk"
include: "rules/pangraph.smk"
include: "rules/resistance.smk"
include: "rules/distances.smk"
include: "rules/assembly_qc.smk"
include: "rules/gubbins.smk"
include: "rules/backbone_joints.smk"
include: "rules/plasmids.smk"
include: "rules/figs.smk"
include: "rules/annotations.smk"
include: "rules/panx.smk"
include: "rules/rates.smk"
include: "rules/hotspots.smk"


rule all:
    input:
        expand(rules.FG_metadata.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_tree_summary.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_homoplasies.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_recombination_filter.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_block_distr_fig.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_distances.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_coresynt.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_circle_synteny.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_junctions_survey.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_junctions_stats.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_junctions_ann.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_rates.output, dset=dset_names, opt=kernel_opts),


localrules:
    download_gbk,
    GM_download_db,
    QC_busco_download,
    Dfinder_models_download,
    PX_download_repo,
