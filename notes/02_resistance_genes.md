# resistance genes

I use [`abricate`](https://github.com/tseemann/abricate) to identify resistance genes in the assemblies.

## conda environment

I installed [`abricate`](https://github.com/tseemann/abricate) following the instructions on the repo:
```bash
conda install -c conda-forge -c bioconda -c defaults abricate
abricate --check
abricate --list
```
This installed `abricate` version `1.0.1`.

## running abricate

`abricate` can be run with different database options:
```
DATABASE	    SEQUENCES	DBTYPE	DATE
ecoli_vf	    2701	nucl	2021-Mar-27
ncbi	        5386	nucl	2021-Mar-27
vfdb	        2597	nucl	2021-Mar-27
megares	        6635	nucl	2021-Mar-27
argannot	    2223	nucl	2021-Mar-27
plasmidfinder	460	    nucl	2021-Mar-27
ecoh	        597	    nucl	2021-Mar-27
resfinder	    3077	nucl	2021-Mar-27
card	        2631	nucl	2021-Mar-27
```

It seems that `ncbi` is more focused on genes that confer resistance by their presence alone, while `card` is more comprehensive and include genes that are somewhat involved in resistance pathways (e.g. efflux pumps...).
I run it with both databases to have both a focused and broad view of the resistance genes.
```bash
abricate -db card isolate.gbk > card/isolate.tab
abricate -db ncbi isolate.gbk > card/isolate.tab
```

I then combine the results with:
```bash
abricate --summary card/*.tab > card/summary.txt
abricate --summary ncbi/*.tab > ncbi/summary.txt
```

