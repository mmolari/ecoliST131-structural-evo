
ref_genome.fa:
	@echo "link file $@"
	ln -s ../../data/fa/NZ_CP124372.1.fa $@

ref_genome.fa.fai: ref_genome.fa
	@echo "index reference genome"
	samtools faidx $<

long_reads.fa:
	@echo "extract long reads"
	seqkit seq reads.fa --min-len 5000 > long_reads.fa

aln.sam: ref_genome.fa long_reads.fa
	@echo "mapping reads to reference genome"
	minimap2 --secondary=no -ax map-ont $< $< > $@

aln.bam: aln.sam
	@echo "convert sam to bam"
	samtools view -b $< > $@

aln_sorted.bam: aln.bam
	@echo "sort bam file"
	samtools sort $< > $@

aln_sorted.bam.bai: aln_sorted.bam
	@echo "index bam file"
	samtools index $<

all: aln_sorted.bam.bai ref_genome.fa.fai

clean:
	rm -f ref_genome.fa long_reads.fa aln.sam aln.bam aln_sorted.bam aln_sorted.bam.bai ref_genome.fa.fai