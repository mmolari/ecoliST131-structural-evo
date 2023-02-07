configfile: "config/config.yml"


include: "rules/downloads.smk"
include: "rules/pangraph.smk"
include: "rules/distances.smk"


localrules:
    download_gbk,
    PG_all,
