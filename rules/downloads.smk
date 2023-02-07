rule download_gbk:
    output:
        "data/gbk/{acc}.gbk",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        ncbi-acc-download {wildcards.acc} -e all -F genbank
        mv {wildcards.acc}.gbk {output}
        """


rule gbk_to_fa:
    input:
        gbk=rules.download_gbk.output,
    output:
        fa="data/fa/{acc}.fa",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/utils/gbk_to_fa.py --gbk {input.gbk} --fa {output.fa}
        """


# datasets_fnames = config["datasets"]

# # read accession numbers from dataset files
# datasets = {}
# for k, fname in datasets_fnames.items():
#     with open(fname, "r") as f:
#         acc_nums = f.readlines()
#     acc_nums = [an.strip() for an in acc_nums]
#     acc_nums = [an for an in acc_nums if len(an) > 0]
#     datasets[k] = acc_nums
# rule download_all:
#     input:
#         expand(rules.gbk_to_fa.output, acc=datasets["ST131"]),
