# Split the big input file into lines of 100, removing empty lines
# could do this with python...
from math import floor
import subprocess
import glob
import os
from time import sleep

APIKEY = "abf81216999e19fc8a23ae76fb9c6b208908" # my NCBI API key

id_dict = {}
with open('id_and_biosample.csv', 'r') as f:
	for i, line in enumerate(f.readlines()):
		nuccore_id, biosample_id = line.strip().split(',')
		slot = floor(i/100)
		if slot in id_dict.keys():
			id_dict[slot][nuccore_id] = biosample_id
		else:
			id_dict[slot] = {nuccore_id : biosample_id }

#for i in range(len(id_dict.keys())):
for i in range(len(id_dict.keys())):
	tmp_file = 'tmp_biosamples_'+str(i)+'.txt'
	with open(tmp_file, 'w') as f:
		for nuccore_id, biosample_id in id_dict[i].items():
			f.write('%s\n' % biosample_id)
	with open(tmp_file+'.xml', 'w') as f:
		print('reading in tmp_file:', tmp_file)
		epost_command = ['epost', '-db', 'biosample', '-input', tmp_file, 
				'api_key', APIKEY]	
		epost_process = subprocess.run(epost_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
		docsum_output = subprocess.run(['efetch', '-format', 'docsum'],
                              input=epost_process.stdout, stdout=f, stderr=subprocess.PIPE)	
		print('executed epost command')
	#sleep(10) 


# Download with epost
# Run the python extract-xml.py

