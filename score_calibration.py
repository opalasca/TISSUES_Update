#!/usr/bin/python
import sys,getopt,re, gzip
import subprocess
import math
from collections import defaultdict
import matplotlib.pyplot as plt
plt.rcParams.update({'axes.labelsize': 'large'})
plt.rcParams.update({'axes.titlesize': 'x-large'})
from os import walk
import gzip
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

datasets=['human_text_mining','mouse_gnf','mouse_gnfv3','mouse_rnaseq_encode','mouse_rnaseq_mit', 'rat_array', 'rat_rnaseq_mit', 'rat_rnaseq_bodymap', 'pig_array', 'pig_rnaseq_aarhus', 'pig_rnaseq_wur']
colours = ["#556b2f","#4169e1","#4682b4","#1e90ff","#191970","#cd853f","#daa520","#a0522d","#ffc0cb","#ffa07a","#db7093"]
names = ["Text mining","Mouse GNF", "Mouse GNF V3", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat Array","Rat RNA-seq MIT","Rat RNA-seq BodyMap", "Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR"]

uniprot_fe = {}
uniprot_score = {}

for i in range(11):
	uniprot_fe[datasets[i]] = []
	uniprot_score[datasets[i]] = []

for i in range(11):
	u = open("data/"+datasets[i]+"_uniprot_orth_fold_enrichment_analysis.tsv")
	print datasets[i]
	for line in u:
		cols = re.split('\t', line)
		uniprot_score[datasets[i]].append(float(cols[0]))
		uniprot_fe[datasets[i]].append(float(cols[1]))
	u.close()

ticksy = np.array([0, 1,2,3])		

fig = plt.figure(figsize=(8,7))
ax=plt.subplot(111)
for i in range(11):
	plt.scatter(uniprot_score[datasets[i]], uniprot_fe[datasets[i]], color=colours[i],label=names[i],s=5)
plt.xlim(-2, 4)
plt.yticks(ticksy)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.legend(loc='center left', markerscale=2., scatterpoints=1, fontsize=11, bbox_to_anchor=(1, 0.5))
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=3, markerscale=2., scatterpoints=1, fontsize=11)
ax.set_axis_bgcolor('none')
ax.set_ylabel('Fold enrichment compared to gold standard')
ax.set_xlabel('Mean confidence score')

plt.show()
fig.savefig('figures/score_calibration.pdf', bbox_inches='tight')
fig.savefig('figures/score_calibration.png', bbox_inches='tight',dpi=100)
plt.close(fig)


