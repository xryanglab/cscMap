# !/usr/bin/python2.7
# -*- coding: utf-8 -*-
# Yuting Wang
# 20210623

'''
The pipeline "cscMap" is identified the cross-strand chimeric RNA(cscRNA) from RNA-seq data.
usage: python2 $Dir/cscRNA_identify.py readW.sam readC.sam r1_need_r2_90M.txt r2_need_r1_90M.txt junctionSite.txt

# File guidance
# r1_need_r2_90M.txt: read1 cross the junction loci
# r2_need_r1_90M.txt: read2 cross the junction loci
# junctionSite.txt: the sum of junction loci

# Requirement
# correct files : refs/Database/hg19.chrom.sizes
#                 refs/hg19.fa.out.bed
# replace the Dir file "../refs/" with your own directory file!!!!
# The pipeline is under the python2.7 environment and required the following softwares:tophat2, bowtie2, samtools, bedtools(intersect), RSEQC(bam2wig.py), BEDOPS(sam2bed).

'''



import re
import os
from sys import argv,stderr
fo1 = argv[1]
fo2 = argv[2]

## getMapReadsForward()
os.system("awk '$2==0' "+fo1+" >readW_map.sam")

## getMapReadsReverse()
os.system("awk '$2==16' "+fo2+"> readC_map.sam")

## remove repeat sequences

def filter_repeats(in1,in2):
	read_noRepeats_r1 = []
	read_noRepeats_r2 = []
	arrayC_set = set()
	with open(in1) as f:
		for line in f:
			tmp = line.split("\t")
			arrayC_set.add("%s" % (tmp[0].strip()))
	with open(in2) as f:
		for line in f:
			tmpW = line.split("\t")
			name = "%s" % (tmpW[0].strip())
			if name in arrayC_set:
				continue
			else:
				if "_2" in name:
					read_noRepeats_r2.append(line)
				else:
					read_noRepeats_r1.append(line)
	return read_noRepeats_r1,read_noRepeats_r2

os.system("sam2bed < readW_map.sam > readW_map.bed")
os.system("cut -f1-6 readW_map.bed >readW_map.bed6")
os.system("bedtools intersect -a readW_map.bed6 -b /Share2/home/yangxr/LabMem/ytwang/Database/hg19.fa.out.bed  -wa -wb > readW_repeats_overlap.txt")
os.system("cut -f 4 readW_repeats_overlap.txt  |sort -u > readW_repeats_overlap_lable.txt")
[readW_noRepeats_r1,readW_noRepeats_r2] = filter_repeats("readW_repeats_overlap_lable.txt","readW_map.sam")

os.system("sam2bed < readC_map.sam > readC_map.bed")
os.system("cut -f1-6 readC_map.bed > readC_map.bed6")
os.system("bedtools intersect -a readC_map.bed6 -b /Share2/home/yangxr/LabMem/ytwang/Database/hg19.fa.out.bed  -wa -wb > readC_repeats_overlap.txt")
os.system("cut -f 4 readC_repeats_overlap.txt  |sort -u > readC_repeats_overlap_lable.txt")
[readC_noRepeats_r1,readC_noRepeats_r2] = filter_repeats("readC_repeats_overlap_lable.txt","readC_map.sam")


def filter_chr(array1,array2):
	array = []
	arr1 = {}
	for i in range(len(array1)):
		tmp = array1[i].split("\t")
		name = "%s%s" % (tmp[0],tmp[2])
		arr1[name]=array1[i]
	for j in range(len(array2)):
		tmp = array2[j].split("\t")
		name2 = "%s%s" % (tmp[0],tmp[2])
		if name2 in arr1:
			array.append(arr1[name2])
			array.append(array2[j])
	return array

array_sameChr_r1 = filter_chr(readW_noRepeats_r1,readC_noRepeats_r1)
array_sameChr_r2 = filter_chr(readW_noRepeats_r2,readC_noRepeats_r2)

def filter_arms(array):
	arr = []
	tmp1 = []
	tmp2 = []
	line1 = ""
	for i,l in enumerate(array):
		flag = l.split("\t")[5]
		nums = re.split('[A-Z]',flag)[:-1]
		chars = re.split('\d+',flag)[1:]
		if i % 2 == 0:
			line1 = l
		for k,c in enumerate(chars):
			if c == 'S' or c == 'M':
				if int(nums[k]) > 20:
					if i % 2 == 0:
						tmp1.append(c)
					else:
						tmp2.append(c)
		if i % 2 == 1:
			s1 = "".join(tmp1)
			s2 = "".join(tmp2)
			if len(s1) > len(s2):
				s1,s2 = s2,s1
			outflag = 0
			for m in range(0,len(s1),2):
				if m == len(s1)-1: break
				if s1[m:(m+2)] in s2 and s1[m] !=s1[m+1]:
					outflag = 1
			if outflag == 1:
				arr.append(line1)
				arr.append(l)
			elif len(tmp1) != len(tmp2):
				stderr.write("warning: the pair is:%s\t%s\n" % (line1.split("\t")[5],l.split("\t")[5]))
			tmp1 = []
			tmp2 = []
			line1 = ""
	return arr

array_r1_need = filter_arms(array_sameChr_r1)
array_r2_need = filter_arms(array_sameChr_r2)

def filter_90M(array):
	arr_dic = {}
	for i in range(len(array)):
		tmpC = array[i].split("\t")
		name = "%s" % (tmpC[0])
		strC=tmpC[5].split("M")
		b=filter(str.isdigit,strC[0][-3:])
		match=int(b)
		if(match>=90):
			arr_dic[name]=array[i]
	return arr_dic

readW_noRepeats_r1.extend(readC_noRepeats_r1)
readW_noRepeats_r2.extend(readC_noRepeats_r2)
array_r1_90M_dic = {}
array_r2_90M_dic = {}
array_r1_90M_dic = filter_90M(readW_noRepeats_r1)
array_r2_90M_dic = filter_90M(readW_noRepeats_r2)

def filter_overlap_reviseLoci(array1,dic2):
	array = []
	for i in range(0,len(array1),2):
		tmp1=array1[i].split("\t")
		if "_2" not in tmp1[0]: # r1_need
			name = tmp1[0] + "_2"
			if name in dic2:
				tmp2 = array1[i+1].split("\t")
				tmp3 = dic2[name].split("\t")
				tmp1_chars = re.split('\d+',tmp1[5])[1:]
				tmp2_chars = re.split('\d+',tmp2[5])[1:]
				tmp1_num = int(filter(str.isdigit,re.split('M',tmp1[5])[0][-2:]))
				tmp2_num = int(filter(str.isdigit,re.split('M',tmp2[5])[0][-2:]))
				tmp3_num = int(filter(str.isdigit,tmp3[5].split("M")[0][-3:]))
				tmp1_s = "".join(tmp1_chars)
				tmp2_s = "".join(tmp2_chars)
				loci1 = int(tmp1[3])
				loci2 = int(tmp2[3])
				loci3 = int(tmp3[3])
				loci1_end = loci1 + tmp1_num
				loci2_end = loci2 + tmp2_num
				loci3_end = loci3 + tmp3_num
				chr2 = tmp2[2]
				chr3 = tmp3[2]
				mismatch1 = int(array1[i][(array1[i].find("NM:i:")+5):(array1[i].find("NM:i:")+7)])
				mismatch2 = int(array1[i+1][(array1[i+1].find("NM:i:")+5):(array1[i+1].find("NM:i:")+7)])
				mismatch3 = int(dic2[name][(dic2[name].find("NM:i:")+5):(dic2[name].find("NM:i:")+7)])
				if chr2 == chr3:
					if mismatch1<3 and mismatch2<3 and mismatch3<3:
						if ("MS" in tmp1_s) and ("MS" in tmp2_s) and (int(tmp1_num)+int(tmp2_num) >=101):
							if(tmp1[1]=="0")and(tmp2[1]=="16")and(tmp3[1]=="16"):
								if(loci3_end < loci1_end):
									array.append(array1[i])
									array.append(array1[i+1])
									array.append(dic2[name])
							elif(tmp1[1]=="0")and(tmp2[1]=="16")and(tmp3[1]=="0"):
								if(loci3_end < loci2_end):
									array.append(array1[i])
									array.append(array1[i+1])
									array.append(dic2[name])
						elif ("SM" in tmp1_s) and ("SM" in tmp2_s) and (int(tmp1_num)+int(tmp2_num) >=101):
							if(tmp1[1]=="0")and(tmp2[1]=="16")and(tmp3[1]=="0"):
								if(loci3 > loci2):
									array.append(array1[i])
									array.append(array1[i+1])
									array.append(dic2[name])
							elif(tmp1[1]=="0")and(tmp2[1]=="16")and(tmp3[1]=="16"):
								if(loci3 > loci1):
									array.append(array1[i])
									array.append(array1[i+1])
									array.append(dic2[name])
		else:
			name = tmp1[0][:-2]
			if name in dic2: # r2_need
				tmp1 = dic2[name].split("\t")
				tmp2 = array1[i].split("\t")
				tmp3 = array1[i+1].split("\t")
				tmp2_chars = re.split('\d+',tmp2[5])[1:]
				tmp3_chars = re.split('\d+',tmp3[5])[1:]
				tmp1_num = int(filter(str.isdigit,tmp1[5].split("M")[0][-3:]))
				tmp2_num = int(filter(str.isdigit,re.split('M',tmp2[5])[0][-2:]))
				tmp3_num = int(filter(str.isdigit,re.split('M',tmp3[5])[0][-2:]))
				tmp2_s = "".join(tmp2_chars)
				tmp3_s = "".join(tmp3_chars)
				loci1 = int(tmp1[3])
				loci2 = int(tmp2[3])
				loci3 = int(tmp3[3])
				loci1_end = loci1 + tmp1_num
				loci2_end = loci2 + tmp2_num
				loci3_end = loci3 + tmp3_num
				chr1 = tmp1[2]
				chr2 = tmp2[2]
				mismatch1 = int(array1[i][(array1[i].find("NM:i:")+5):(array1[i].find("NM:i:")+7)])
				mismatch2 = int(array1[i+1][(array1[i+1].find("NM:i:")+5):(array1[i+1].find("NM:i:")+7)])
				mismatch3 = int(dic2[name][(dic2[name].find("NM:i:")+5):(dic2[name].find("NM:i:")+7)])
				if chr1 == chr2:
					if mismatch1<3 and mismatch2<3 and mismatch3<3:
						if ("MS" in tmp2_s) and ("MS" in tmp3_s) and (int(tmp2_num)+int(tmp3_num) >=101):
							if(tmp1[1]=="16")and(tmp2[1]=="0")and(tmp3[1]=="16"):
								if(loci1_end < loci2_end):
									array.append(array1[i])
									array.append(array1[i+1])
									array.append(dic2[name])
							elif(tmp1[1]=="0")and(tmp2[1]=="0")and(tmp3[1]=="16"):
								if(loci1_end < loci3_end):
									array.append(array1[i])
									array.append(array1[i+1])
									array.append(dic2[name])
						elif ("SM" in tmp2_s) and ("SM" in tmp3_s) and (int(tmp2_num)+int(tmp3_num) >=101):
							if(tmp1[1]=="0")and(tmp2[1]=="0")and(tmp3[1]=="16"):
								if(loci1 > loci3):
									array.append(array1[i])
									array.append(array1[i+1])
									array.append(dic2[name])
							elif(tmp1[1]=="16")and(tmp2[1]=="0")and(tmp3[1]=="16"):
								if(loci1 > loci2):
									array.append(array1[i])
									array.append(array1[i+1])
									array.append(dic2[name])
	return array

array_r1_need_r2_90M = filter_overlap_reviseLoci(array_r1_need,array_r2_90M_dic)
array_r2_need_r1_90M = filter_overlap_reviseLoci(array_r2_need,array_r1_90M_dic)

## for output files
fout1 = argv[3]
fout2 = argv[4]
output1 = open(fout1,"w")
output2 = open(fout2,"w")
for i in range(len(array_r1_need_r2_90M)):
	output1.write(array_r1_need_r2_90M[i])
for i in range(len(array_r2_need_r1_90M)):
	output2.write(array_r2_need_r1_90M[i])
output1.close()
output2.close()


# # prepare files for IGV 
file1 = fout1.split(".")[0]
file2 = fout2.split(".")[0]
file_names=[file1,file2]

for i in range(0,2):
	os.system("head -n27 "+fo1+" |cat - "+file_names[i]+".txt >"+file_names[i]+".sam")
	os.system("samtools view -bS "+file_names[i]+".sam >"+file_names[i]+".bam")
	os.system("samtools sort "+file_names[i]+".bam >"+file_names[i]+"_sort.bam")
	os.system("samtools index "+file_names[i]+"_sort.bam")
	os.system("bam2wig.py -i "+file_names[i]+"_sort.bam -s /Share2/home/yangxr/LabMem/ytwang/Database/hg19.chrom.sizes -o "+file_names[i]+"_sort")

# statistic junction loci and counts
fout3 = argv[5]
output = open(fout3,"w")
loci1_dic={}
loci2_dic={}
for i in range(0,2):
	name=file_names[i]+".txt"
	fo = open(name)
	while True:
		line1 = fo.readline()
		line2 = fo.readline()
		line3 = fo.readline()
		if line1=="":
			break
		tmp1 = line1.split("\t")
		tmp2 = line2.split("\t")
		tmp1_chars = re.split('\d+',tmp1[5])[1:]
		tmp2_chars = re.split('\d+',tmp2[5])[1:]
		tmp1_num = filter(str.isdigit,re.split('M',tmp1[5])[0][-2:])
		tmp2_num = filter(str.isdigit,re.split('M',tmp2[5])[0][-2:])
		tmp1_s = "".join(tmp1_chars)
		tmp2_s = "".join(tmp2_chars)
		loci1 = int(tmp1[3])
		loci2 = int(tmp2[3])
		if ("MS" in tmp1_s) and ("MS" in tmp2_s):
			junctionLoci1 = loci1+int(tmp1_num)-1
			junctionLoci2 = loci2+int(tmp2_num)-1
			if junctionLoci1 > junctionLoci2:
				junctionLoci1,junctionLoci2 = junctionLoci2,junctionLoci1
			for i in range(5):
				lociP = str(int(junctionLoci1)+i+1)
				lociM = str(int(junctionLoci1)-i-1)
				lociP2 = str(int(junctionLoci2)+i+1)
				lociM2 = str(int(junctionLoci2)-i-1)
				str1=tmp1[2]+"\t"+lociP
				str2=tmp1[2]+"\t"+lociM
				str1_2=tmp1[2]+"\t"+lociP2
				str2_2=tmp1[2]+"\t"+lociM2
				loci1_dic[str1]=str(junctionLoci1)
				loci1_dic[str2]=str(junctionLoci1)
				loci2_dic[str1_2]=str(junctionLoci2)
				loci2_dic[str2_2]=str(junctionLoci2)
			output.write("MS:"+"\t"+tmp1[2]+"\t"+str(junctionLoci1)+"\t"+tmp2[2]+"\t"+str(junctionLoci2)+"\n")
		elif ("SM" in tmp1_s) and ("SM" in tmp2_s):
			junctionLoci1=loci1
			junctionLoci2=loci2
			if junctionLoci1 > junctionLoci2:
				junctionLoci1,junctionLoci2 = junctionLoci2,junctionLoci1
			for i in range(5):
				lociP = str(int(junctionLoci1)+i+1)
				lociM = str(int(junctionLoci1)-i-1)
				lociP2 = str(int(junctionLoci2)+i+1)
				lociM2 = str(int(junctionLoci2)-i-1)
				str1=tmp1[2]+"\t"+lociP
				str2=tmp1[2]+"\t"+lociM
				str1_2=tmp1[2]+"\t"+lociP2
				str2_2=tmp1[2]+"\t"+lociM2
				loci1_dic[str1]=str(junctionLoci1)
				loci1_dic[str2]=str(junctionLoci1)
				loci2_dic[str1_2]=str(junctionLoci2)
				loci2_dic[str2_2]=str(junctionLoci2)
			output.write("SM:"+"\t"+tmp1[2]+"\t"+str(junctionLoci1)+"\t"+tmp2[2]+"\t"+str(junctionLoci2)+"\n")

output.close()
# os.system('''sort '''+fout3+''' |uniq -c |sort -r -n |awk -F'[MS|SM]' '$1>=2' |awk '{print $3"\t"$4"\t"$5"\t"$6"\t"$2"\t"$1;}'|sed 's/://' > junctionLibrary_count.txt''')

# delete the tmp files
os.system("rm readC_map.*")
os.system("rm readW_map.*")
os.system("rm readC_repeats_overlap*")
os.system("rm readW_repeats_overlap*")
os.system("rm r1_need_r2_90M.sam")
os.system("rm r1_need_r2_90M.bam")
os.system("rm r2_need_r1_90M.sam")
os.system("rm r2_need_r1_90M.bam")

