# ST43 / ST131

## dataset preparation

The accession numbers were obtained by Liam with the following procedure:

1. downloaded all fasta (02.02.2023) from the following search in RefSeq: `st131 escherichia coli AND srcdb_refseq[PROP]`. With length filters 4-7Mbp to remove any bad genomes. This gives n=478 genome hits, not all of which contain sequencing data (because people embargoed the data).Downloading fasta we have n=133 genomes in single fasta.
2. Splitting fasta and ignoring duplicate genome IDs there were n=109 genomes. Liam applied [mlst](https://github.com/tseemann/mlst) to check the sequence types (ST)
```bash
for f in *fa;
do
	mlst --quiet --scheme ecoli $f >> mlst-results.tsv
	echo $f
done 
```
ST131 is called ST43 in the standard E.Coli MLST scheme. Filtering for these:
```bash
awk '$3==43' mlst-results.tsv | awk '{print $1}' | cut -d '.' -f 1 > ST43.txt
```
Gives n=86 chromosome accessions which are definitely ST43/ST131. Other minor STs (or untyped) may well also be within the ST43/ST131 group (including even more closely related) but just have differences at the alleles that define the ST. This might be a good dataset with expected evolutionary divergence <100 years to investigate.