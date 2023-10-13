rule QC_busco_download:
    output:
        directory("data/busco_models"),
    conda:
        "../conda_env/busco.yml"
    shell:
        """
        busco --download enterobacterales_odb10
        mv busco_downloads {output}
        """


# busco --download prokaryota --download_path {output}


rule QC_busco_run:
    input:
        fa=rules.gbk_to_fa.output.fa,
        mod=rules.QC_busco_download.output,
    output:
        directory("data/busco/{acc}"),
    conda:
        "../conda_env/busco.yml"
    shell:
        """
        busco -i {input.fa} \
            --offline \
            -l enterobacterales_odb10 \
            --download_path {input.mod} \
            -m genome \
            -o {output}
        """


rule QC_mlst:
    input:
        gbks=lambda w: expand(
            rules.download_gbk.output,
            acc=sorted(dset_chrom_accnums[w.dset] + excluded[w.dset]),
        ),
    output:
        "results/{dset}/assembly_qc/mlst/{scheme}.tsv",
    conda:
        "../conda_env/mlst.yml"
    shell:
        """
        mlst {input.gbks} \
            --scheme {wildcards.scheme} \
            > {output}
        """


rule QC_summary:
    input:
        lambda w: expand(
            rules.QC_busco_run.output,
            acc=sorted(dset_chrom_accnums[w.dset] + excluded[w.dset]),
        ),
    output:
        "results/{dset}/assembly_qc/busco_summary.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/assembly_qc/busco_summary.py \
            --fld {input} \
            --csv {output}
        """


rule QC_alleles_db:
    input:
        fa=lambda w: config["allele-files"][w.allele],
    output:
        db=directory(
            "data/blast_alleles_db/{allele}",
        ),
    conda:
        "../conda_env/mapping.yml"
    shell:
        """
        gzip -c -d {input.fa} | \
        makeblastdb \
            -title {wildcards.allele} \
            -hash_index \
            -dbtype nucl \
            -parse_seqids \
            -out {output.db}/{wildcards.allele}
        """


rule QC_alleles_map:
    input:
        fa=rules.gbk_to_fa.output.fa,
        db=rules.QC_alleles_db.output.db,
    output:
        "data/alleles/map/{allele}/{acc}.tsv",
    params:
        min_id=95,
    conda:
        "../conda_env/mapping.yml"
    shell:
        """
        blastn -db {input.db}/{wildcards.allele} \
        -outfmt '6 sseqid slen length nident qseqid qstart qend sstrand' \
            -ungapped -dust no -word_size 32 -max_target_seqs 10000 \
            -perc_identity {params.min_id} -evalue 1E-20 \
            -query {input.fa} \
           > {output}
        """
        # -outfmt '6 sseqid slen length nident qseqid qstart qend qseq sstrand' \


rule QC_alleles_assign:
    input:
        paf=rules.QC_alleles_map.output,
    output:
        "data/alleles/assign/{allele}/{acc}.tsv",
    params:
        min_cov=0.99,
        min_id=0.99,
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/assembly_qc/assign_alleles.py \
            --blast_tsv {input.paf} \
            --min_cov {params.min_cov} \
            --min_id {params.min_id} \
            --tsv_out {output}
        """


rule QC_alleles_concat:
    input:
        lambda w: expand(
            rules.QC_alleles_assign.output,
            allele=w.allele,
            acc=sorted(dset_chrom_accnums[w.dset] + excluded[w.dset]),
        ),
    output:
        "results/{dset}/assembly_qc/alleles/{allele}.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        echo "iso\tlocus\tmatch\tallele\tcoverage\tsim\tmatches\taln_L\tallele_L" > {output}
        cat {input} >> {output}
        """


rule QC_alleles_summary:
    input:
        dfs=expand(
            rules.QC_alleles_concat.output,
            allele=config["chr_alleles"],
            allow_missing=True,
        ),
    output:
        csv="results/{dset}/assembly_qc/alleles_summary.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/assembly_qc/alleles_summary.py \
            --in_dfs {input.dfs} \
            --summary_df {output.csv}
        """


rule QC_all:
    input:
        expand(rules.QC_summary.output, dset=dset_names),
        expand(rules.QC_mlst.output, dset=dset_names, scheme=config["mlst_schemes"]),
        expand(rules.QC_alleles_summary.output, dset=dset_names),
