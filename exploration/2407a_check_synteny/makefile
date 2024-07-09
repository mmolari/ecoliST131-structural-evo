
ACC = \
	NZ_CP049085.2 \
	NZ_CP124320.1 \
	NZ_CP124328.1 \
	NZ_CP124372.1 \
	NZ_CP124481.1 \
	NZ_OW968277.1 \
	NZ_OX637959.1 \
	NZ_OX637961.1

NEIGH = \
	NZ_CP049085.2__NZ_CP103557.1 \
	NZ_CP124320.1__NZ_CP124500.1 \
	NZ_CP124328.1__NZ_CP124355.1 \
	NZ_CP124372.1__NZ_CP107137.1 \
	NZ_CP124481.1__NZ_CP124487.1 \
	NZ_OW968277.1__NZ_CP107147.1 \
	NZ_OX637959.1__NZ_CP103742.1 \
	NZ_OX637961.1__NZ_CP103742.1

# Function to get value from the lookup table
find_nb = $(word 2, $(subst __, , $(filter $1__%,$(NEIGH))))

genomes:
	mkdir -p genomes

genomes/%.fa: genomes
	ln -sf ../../../data/fa/$*.fa $@

genomes/%.fa.fai: genomes/%.fa
	samtools faidx $<

nb_aln:
	mkdir -p nb_aln

nb_aln/%.paf: nb_aln genomes/%.fa
	value=$(call find_nb,$*); \
	$(MAKE) "genomes/$${value}.fa"; \
	minimap2 -x asm5 "genomes/$${value}.fa" genomes/$*.fa > $@

read_aln:
	mkdir -p read_aln

read_aln/%.sam: genomes/%.fa reads/%.fa.gz read_aln
	seqkit seq reads/$*.fa.gz --min-len 5000 > long_reads_$*.fa
	minimap2 -ax map-ont genomes/$*.fa long_reads_$*.fa > $@
	rm long_reads_$*.fa

read_aln/%.bam: read_aln/%.sam
	samtools view -b $< > $@

read_aln/%_sorted.bam: read_aln/%.bam
	samtools sort $< > $@

read_aln/%_sorted.bam.bai: read_aln/%_sorted.bam
	samtools index $<

.PRECIOUS: \
	genomes/%.fa \
	genomes/%.fa.fai \
	nb_aln/%.paf \
	read_aln/%_sorted.bam \
	read_aln/%_sorted.bam.bai

all:
	$(MAKE) $(addprefix genomes/, $(addsuffix .fa.fai, $(ACC)))
	# $(MAKE) $(addprefix nb_aln/, $(addsuffix .paf, $(ACC)))
	# $(MAKE) $(addprefix read_aln/, $(addsuffix _sorted.bam.bai, $(ACC)))

clean:
	rm -rf genomes nb_aln read_aln