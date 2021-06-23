# !/usr/bin/python2 by ytwang in 2016/11/7
# This script remove the chimeric RNA reads mapped to the transcript.
# Usage: python rmToTrspt.py r1_need_r2_90M.txt ../tophat_normalMapping1/unmapped.lables r1_need_r2_90M_f.txt

from sys import argv,stderr
inF1 = argv[1]
inF2 = argv[2]
fout = argv[3]
output = open(fout,"w")
unmap_set=set()
with open(inF2) as f:
	for line in f:
		tmp = line.split("\t")
		unmap_set.add("%s" % (tmp[0].strip()))
		
	

fo = open(inF1)
while True:
	line1 = fo.readline()
	line2 = fo.readline()
	line3 = fo.readline()
	if line1=="":
		break
	tmp = line1.split("\t")[0].strip()
	if (tmp in unmap_set):
		output.write(line1)
		output.write(line2)
		output.write(line3)
		
	

fo.close()
output.close()
