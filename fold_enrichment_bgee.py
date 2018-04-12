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

datasets=['human_all_data_bgee_ranking','human_integrated','mouse_all_data_bgee_ranking','mouse_integrated','rat_all_data_bgee_ranking','rat_integrated', 'pig_all_data_bgee_ranking','pig_integrated']
colours = ["#4169e1", "#4682b4", "#cd853f","#ffa07a","#db7093","#daa520","#4682b4","#cd853f"]
names = ["Human BGEE","Human TISSUES" ,"Mouse BGEE", "Mouse TISSUES","Rat BGEE","Rat TISSUES","Pig BGEE","Pig TISSUES"] 

nb=8

uniprot_fe = {}
uniprot_expr = {}

for i in range(nb):
	uniprot_fe[datasets[i]] = []
	uniprot_expr[datasets[i]] = []

for i in range(nb):
	u = open("data/"+datasets[i]+"_uniprot_orth_fold_enrichment_analysis.tsv")
	#u = open("data/"+datasets[i]+"_uniprot_fold_enrichment_analysis.tsv")
	print datasets[i]
	for line in u:
		cols = re.split('\t', line)
		if (float(cols[2])==0):
			expr=0.0001
		else:
			expr=float(cols[2])	
		uniprot_expr[datasets[i]].append(expr)
		uniprot_fe[datasets[i]].append(float(cols[1]))
	u.close()

def newplot(ax, data_expr, data_fe, col, name, low, medium, high):
	plt.scatter(data_expr, data_fe, color = col, label=name)
	#plt.xscale('log')
	#ax.get_xaxis().set_visible(False)
#	plt.ylim((0, limity))
#	plt.xlim(0.0005, 300000)
	plt.yticks(ticksy)
	ax.text(textx, texty, name, color=col, verticalalignment='center', horizontalalignment='left',transform=ax.transAxes, fontsize=12)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	#ax.xaxis.set_ticks_position('bottom')

left = -3
right = 4
limity = 4
ticksy = np.array([0, 1,2,3])		
textx=0.03
texty=0.8


fig = plt.figure(figsize=(9,12))

for i in range(nb):
	k=i+1;
	ax = plt.subplot(nb,1,k)
	newplot(ax, uniprot_expr[datasets[i]], uniprot_fe[datasets[i]], colours[i],names[i], 40, 50,100)
ax.get_xaxis().set_visible(True)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

big_ax = fig.add_subplot(111)
big_ax.set_axis_bgcolor('none')
big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
big_ax.set_ylabel('Fold enrichment')
big_ax.set_xlabel('Mean expression rank score')

plt.show()
fig.savefig('figures/fold_enrichment_bgee_unfiltered.pdf', bbox_inches='tight')
fig.savefig('figures/fold_enrichment_bgee_unfiltered.png', bbox_inches='tight',dpi=100)
plt.close(fig)


