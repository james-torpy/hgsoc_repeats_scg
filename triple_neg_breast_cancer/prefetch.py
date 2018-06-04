### prefetch.py ###
# downloads sra files and calls Snakefile:

import os
import re


#############################################################
### 0. Set up variables and paths ###
#############################################################

exp_name = 'SAMN03979292'


home_dir = '/share/ScratchGeneral/jamtor/'

project_dir = home_dir + \
	'/projects/hgsoc_repeats/RNA-seq/triple_neg_breast_cancer/'

raw_dir = project_dir + '/raw_files/' + exp_name

results_dir = project_dir + '/results/' + exp_name

sra_dir = raw_dir + '/sra/'

if not os.path.exists(project_dir + '/' + sra_dir):
	os.makedirs(project_dir + '/' + sra_dir)


#############################################################
### 1. Fetch input variables and prefetch files ###
#############################################################

# fetch sra accessions:
with open(raw_dir + '/SraAccList.txt') as f:
    accessions = f.read().splitlines()

# remove empty strings from accessions list:
SAMPLES = list(filter(None, accessions))

#SAMPLES = ['SRR600983-input', 'SRR600956-H3K4me3', 
#	'SRR600559-H3K27me3']

print('SAMPLES are: ' + ' '.join(SAMPLES))

# prefetch sample sras if they don't already exist:
for i, s in enumerate(SAMPLES):
	if not os.path.isfile(sra_dir + s + ".sra"):
		print('Prefetching ' + s)
		os.system('prefetch ' + s)
		os.system('mv /home/jamtor/ncbi/public/sra/' + s + \
			'.sra ' + project_dir + sra_dir)
	else:
		print(s + ' already prefetched')

## call Snakemake
os.system('snakemake -s ' + project_dir + \
	'repeat_quantification_snakefile')
