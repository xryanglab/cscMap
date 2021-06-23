# !/usr/bin/python2.7 by ytwang in 2016/11/7
# This script get the chimeric RNA fastq
# Usage: python ../getFQ.py r1_need_r2_90M.txt r1_need_r2_90M.fq

from sys import argv,stderr
inF1 = argv[1]
fout = argv[2]
output = open(fout,"w") 
fo = open(inF1)
while True:
	line1 = fo.readline()
	line2 = fo.readline()
	line3 = fo.readline()
	if line1=="":
		break
	tmp = line1.split("\t")
	output.write("@"+tmp[0]+"\n"+tmp[9]+"\n"+"+"+"\n"+tmp[10]+"\n")
	
