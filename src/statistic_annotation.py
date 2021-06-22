# !/usr/bin/python2.7
# -*- coding: utf-8 -*-
# Yuting Wang
# 20201014

'''
The script is used to annotate the junction sites.
usage: python ../statistic_annotation.py */Database/gencode.v19.annotation_noHead.gtf loci.txt loci1_annotation5bp+gtf.txt loci2_annotation5bp+gtf.txt

'''

#
from sys import argv,stderr
inF1 = argv[1]
inF2 = argv[2]
# fo1 = argv[3]
# fo2 = argv[4]
fout1 = argv[3]
fout2 = argv[4]
# out1 = open(fo1,"w")
# out2 = open(fo2,"w")
output1 = open(fout1,"w")
output2 = open(fout2,"w")
fo = open(inF2)
array = []
with open(inF1) as f:
	for line in f:
		array.append(line.split("\t"))
		
	

while True:
	line = fo.readline()
	if line:
		tmp   = line.split("\t")
		chrA  = tmp[0]
		loci1 = int(tmp[1])
		loci2 = int(tmp[2])
		flag  = tmp[3].strip()
		distance1 = 0
		distance2 = 0
		type1_set = set()
		type2_set = set()
		for i in range(len(array)):
			chrB = array[i][0]
			gene_type = array[i][2]
			star = int(array[i][3])
			end  = int(array[i][4])
			if(chrA==chrB):
				if((loci1<=end)and(loci1>=star))or((loci1+5<=end)and(loci1+5>=star))or((loci1-5<=end)and(loci1-5>=star)):
					type1_set.add(gene_type)
					if gene_type=="exon":
						if flag=="SM":
							distance1 = abs(star-loci1)
							output1.write(chrA+"\t"+str(loci1)+"\t"+array[i][0]+"\t"+array[i][1]+"\t"+array[i][2]+"\t"+array[i][3]+"\t"+array[i][4]+"\t"+array[i][5]+"\t"+array[i][6]+"\t"+array[i][7]+"\t"+array[i][8].strip()+"\t"+"distance:"+str(distance1)+"\n")
						else:
							distance1 = abs(end-loci1)
							output1.write(chrA+"\t"+str(loci1)+"\t"+array[i][0]+"\t"+array[i][1]+"\t"+array[i][2]+"\t"+array[i][3]+"\t"+array[i][4]+"\t"+array[i][5]+"\t"+array[i][6]+"\t"+array[i][7]+"\t"+array[i][8].strip()+"\t"+"distance:"+str(distance1)+"\n")
					else:
						output1.write(chrA+"\t"+str(loci1)+"\t"+array[i][0]+"\t"+array[i][1]+"\t"+array[i][2]+"\t"+array[i][3]+"\t"+array[i][4]+"\t"+array[i][5]+"\t"+array[i][6]+"\t"+array[i][7]+"\t"+array[i][8])
				if((loci2<=end)and(loci2>=star))or((loci2+5<=end)and(loci2+5>=star))or((loci2-5<=end)and(loci2-5>=star)):
					type2_set.add(gene_type)
					if gene_type=="exon":
						if flag=="SM":
							distance2 = abs(star-loci2)
							output2.write(chrA+"\t"+str(loci2)+"\t"+array[i][0]+"\t"+array[i][1]+"\t"+array[i][2]+"\t"+array[i][3]+"\t"+array[i][4]+"\t"+array[i][5]+"\t"+array[i][6]+"\t"+array[i][7]+"\t"+array[i][8].strip()+"\t"+"distance:"+str(distance2)+"\n")
						else:
							distance2 = abs(end-loci2)
							output2.write(chrA+"\t"+str(loci2)+"\t"+array[i][0]+"\t"+array[i][1]+"\t"+array[i][2]+"\t"+array[i][3]+"\t"+array[i][4]+"\t"+array[i][5]+"\t"+array[i][6]+"\t"+array[i][7]+"\t"+array[i][8].strip()+"\t"+"distance:"+str(distance2)+"\n")
					else:
						output2.write(chrA+"\t"+str(loci2)+"\t"+array[i][0]+"\t"+array[i][1]+"\t"+array[i][2]+"\t"+array[i][3]+"\t"+array[i][4]+"\t"+array[i][5]+"\t"+array[i][6]+"\t"+array[i][7]+"\t"+array[i][8])
					if("gene" in type1_set)and("transcript" in type1_set)and("exon" in type1_set)and("gene" in type2_set)and("transcript" in type2_set)and("exon" in type2_set):
						break
		# if("gene" in type1_set)and("transcript" in type1_set)and("exon" in type1_set):
			# out1.write(chrA+"\t"+str(loci1)+"\t"+"exon"+"\n")
		# elif("gene" in type1_set)and("transcript" in type1_set):
			# out1.write(chrA+"\t"+str(loci1)+"\t"+"intron"+"\n")
		# else:
			# out1.write(chrA+"\t"+str(loci1)+"\t"+"intergenic region"+"\n")
		# if("gene" in type2_set)and("transcript" in type2_set)and("exon" in type2_set):
			# out2.write(chrA+"\t"+str(loci2)+"\t"+"exon"+"\n")
		# elif("gene" in type2_set)and("transcript" in type2_set):
			# out2.write(chrA+"\t"+str(loci2)+"\t"+"intron"+"\n")
		# else:
			# out2.write(chrA+"\t"+str(loci2)+"\t"+"intergenic region"+"\n")
	else:
		break

output1.close()
output2.close()
# out1.close()
# out2.close()