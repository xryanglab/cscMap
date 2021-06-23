# cscMap
cscMap is a bioinformatics pipeline to search for the RNA chimeras resulted from fusions of the transcripts encoded by the two opposite DNA strands.
### 1. Python dependency
cscMap was written in Python 2.7.16 and has dependencies for popular python libraries:

* re
* os

Also, it depends on some bioinformatics packages.

* TopHat v2
* Bowtie v2
* samtools
* bedtools(intersect)
* RSEQC(bam2wig.py)
* BEDOPS(sam2bed)

Tophat, bowtie are used to align the reads to genome and transcriptome.
RSEQC(bam2wig.py), BEDOPS(sam2bed), bedtools(intersect) are used to file format conversion.

### 2. References file
The reference genome of human (hg19) and mouse (mm10) were downloaded from the UCSC Genome Browser, and the reference genomes of other species were from the Ensemble Genome Browser (zebrafish: GRCz11, C. elegans: WBcel235, fruit fly: BDGP6, Saccharomyces: R64-1-1, E.coli: Escherichia_coliK-12_substr.MG1655). 

The genome annotation files of both human and mouse were obtained from GENCODE (human: v19, mouse: vM17), and the annotation files of all the other species were from the Ensemble Genome Browser. 

The pipeline is required the extra files, for example hg19.chrom.sizes and hg19.fa.out.bed in human, which can get these files from the UCSC Genome Browser(http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/).
* hg19.chrom.sizes - Two-column tab-separated text file containing assembly sequence names and sizes.
* hg19.fa.out.bed - RepeatMasker .out file. Jan 29 2009 (open-3-2-7) version of RepeatMasker, RepBase library: RELEASE 20090120

### 3. Run
The first step is used the tophat2 aligment the RNA-set data.
```
sed 's/ /_2 /' FASTQ2|cat FASTQ1 - > merge.fastq ### distinguish the paired-end reads
tophat2 -p 8 --library-type fr-firststrand -G refs/gencode.v19.annotation.gtf -o tophat_normalMapping/ refs/Bowtie2Index/genome merge.fastq
After this step we get the unmapped.sam file, which is converted to unmapped.fastq and as the input for the cscMap pipeline.
```
The second step is executed the cscMap pipeline to identify the cscRNAs.
```
cd example/
bowtie2 -p 8 -q unmapped.fastq -x refs/Bowtie2Index/genome --nofw --local -S readC.sam
bowtie2 -p 8 -q unmapped.fastq -x refs/Bowtie2Index/genome --norc --local -S readW.sam
python src/cscRNA_identify.py readW.sam readC.sam OutputFile/r1_need_r2_90M.txt OutputFile/r2_need_r1_90M.txt OutputFile/junctionSite.txt
```
The last step is annotated the junction sites and statistic analysis for cscRNAs.
```
cd OutputFile/
python src/statistic_annotation.py refs/gencode.v19.annotation_noHead.gtf junctionSite.txt loci1_annotation5bp+gtf.txt loci2_annotation5bp+gtf.txt
python src/statistic_annotation_overlap.py loci1_annotation5bp+gtf.txt loci2_annotation5bp+gtf.txt junctionLibrary_count_f.txt loci_annotation5bp+gtf.txt
```
## Copyright and license
Copyright (c) 2020 Yuting Wang
The cscMap codes are licensed under THU.

