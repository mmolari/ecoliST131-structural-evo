SHELL:=/bin/bash

pangraph:
	source $$(conda info --base)/etc/profile.d/conda.sh && \
	conda activate pangraph && \
	pangraph build --circular -k minimap2 -s 20 -a 100 -b 5 -l 100 \
		data/genomes/*.fa > data/pangraph.json
	
