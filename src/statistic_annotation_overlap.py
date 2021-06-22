# !/usr/bin/python2.7
# -*- coding: utf-8 -*-
# Yuting Wang
# 20201014

'''
The script is used to annotate the junction sites.
usage: python2 ../statistic_annotation_overlap.py loci1_annotation5bp+gtf.txt loci2_annotation5bp+gtf.txt junctionSite.txt loci_annotation5bp+gtf.txt

'''



from sys import argv,stderr
inF1 = argv[1]
inF2 = argv[2]
inF3 = argv[3]
outF = argv[4]
output = open(outF,"w")
fo = open(inF3)
loci_gtf = {}
with open(inF1) as f:
	for line in f:
		if "distance" in line:
			tmp = line.split("\t")[1]
			if tmp in loci_gtf:
				dis = line.split("\t")[11].split(":")[1].strip()
				dis_b=loci_gtf[tmp].split("\t")[11].split(":")[1].strip()
				if int(dis) < int(dis_b):
					loci_gtf[tmp]=line
			else:
				loci_gtf[tmp]=line
				
			
		
	

with open(inF2) as f:
	for line in f:
		if "distance" in line:
			tmp = line.split("\t")[1]
			if tmp in loci_gtf:
				dis = line.split("\t")[11].split(":")[1].strip()
				dis_b=loci_gtf[tmp].split("\t")[11].split(":")[1].strip()
				if int(dis) < int(dis_b):
					loci_gtf[tmp]=line
			else:
				loci_gtf[tmp]=line
				
			
		
	

blank=" "+"\t"+" "+"\t"+" "+"\t"+" "+"\t"+" "+"\t"+" "+"\t"+" "+"\t"+" "+"\t"+" "+"\t"+" "
while True:
	line = fo.readline()
	if line:
		tmp   = line.split("\t")
		chr=tmp[0]
		loci1 = tmp[1]
		loci2 = tmp[2]
		if(loci1 in loci_gtf)and(loci2 in loci_gtf):
			output.write(tmp[3].strip()+"\t"+loci_gtf[loci1].strip()+"\t"+loci_gtf[loci2])
		elif(loci1 not in loci_gtf)and(loci2 in loci_gtf):
			output.write(tmp[3].strip()+"\t"+chr+"\t"+loci1+"\t"+blank+"\t"+loci_gtf[loci2])
		elif(loci1 in loci_gtf)and(loci2 not in loci_gtf):
			output.write(tmp[3].strip()+"\t"+loci_gtf[loci1].strip()+"\t"+chr+"\t"+loci2+"\t"+blank+"\n")
		elif(loci1 not in loci_gtf)and(loci2 not in loci_gtf):
			output.write(tmp[3].strip()+"\t"+chr+"\t"+loci1+"\t"+blank+"\t"+chr+"\t"+loci2+"\t"+blank+"\n")
	else:
		break

fo.close()
output.close()
