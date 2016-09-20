######hw1######

# ahcg_pipeline
Variant calling pipeline for genomic data analysis

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq

#1.copy the pipeline
search ahcg_pipeline in github then folk it.
https://github.com/shashidhar22/ahcg_pipeline


#2. Changing a remote's URL
instruction: https://help.github.com/articles/changing-a-remote-s-url/

git remote set-url https://github.com/cguo74/ahcg_pipeline.git
git remote -v
####below means that the Changing a remote's URL wasn’t success####
origin	https://github.com/shashidhar22/ahcg_pipeline (fetch)
origin	https://github.com/shashidhar22/ahcg_pipeline (push)
####we add origin####
git remote set-url origin https://github.com/cguo74/ahcg_pipeline.git
git remote -v
####below means that the Changing a remote's URL was success####
origin	https://github.com/cguo74/ahcg_pipeline.git (fetch)
origin	https://github.com/cguo74/ahcg_pipeline.git (push)


#3 download emacs(easy to use than vim)
sudo apt-get update
sudo apt-get install emacs24


#4 ignore the unnecessary folder, add any folder under current position
basespace@ubuntu:~/ahcg_pipeline$ ls -a 
emacs .gitignore


#5.Fasta index using Samtools faidx

###download samtools will produce hg19.fa.fai###
sudo apt-get install samtools
instructions: open the page: http://www.htslib.org/doc/samtools.html

####final commend####
samtools faidx hg19.fa


#6.Genome dict file using picard will produce hg19.dict
###download java and picard###
sudo apt-get install openjdk-8-jre-headless
sudo apt-get install picard
###also download picard jar here https://github.com/broadinstitute/picard/releases/tag/2.6.0###
java -jar picard.jar
instructions: https://broadinstitute.github.io/picard/command-line overview.html#CreateSequenceDictionary

####final commend####
java -jar picard.jar CreateSequenceDictionary R=hg19.fa O=reference.dict


#7 something need to correct before running.In the folder named wantignore I have all the files:patient genome-test_r1.fastq; references genome after bowtie processed- hg19.1.bt2; 
#In resources/genome, I have hg19.fa which is original file, hg19.fa.fai which generated from samtool, hg19.dict generated from picard
mv reference.dict ../resources/genome
mv hg19.fa.fai ../resources/genome
gunzip ./resources/dbsnp/dbsnp_138.hg19.vcf.gz

#try run the pipeline
#To access help use the following command:
python3 ahcg_pipeline.py -h
#you will see:


#########final commend########
python3 ahcg_pipeline.py -t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b ./lib/bowtie2-2.2.9/bowtie2 -p ./lib/picard.jar -g ./lib/GenomeAnalysisTK.jar -i ./wantignore/test_r1.fastq ./wantignore/test_r2.fastq -w ./wantignore/hg19 -d ./resources/dbsnp/dbsnp_138.hg19.vcf -r ./resources/genome/hg19.fa -a ./lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o ./hw1

AR NM_000044.3 313700 Androgen receptor 10661 8
ATM NM_000051.3 607585 Serine-protein kinase ataxia telangiectasia mutated 13147 68
BARD1 NM_000465.3 601593 BRCA associated RING domain 1 5523 11
BRCA1 NM_007298.3 113705 Breast cancer 1, early onset 3699 22
BRCA2 NM_000059.3 600185 Breast cancer 2, early onset 11386 27
BRIP1 NM_032043.2 605882 BRCA1-interacting protein 8166 20
CASP8 NM_001080124.1 601763 Apoptosis-related cysteine protease 8 2750 9
CDH1 NM_004360.3 192090 Cadherin 1 4815 16
CHEK2 NM_001005735.1 604373 Serine/threonine checkpoint kinase 2 1991 16
DIRAS3 NM_004675.2 605193 GTP-binding Ras-like protein 3 1642 2
ERBB2 NM_001005862.1 164870 Avian erythroblastic leukemia viral oncogene homolog 2 4816 30
NBN NM_002485.4 602667 Nibrin 4639 16
PALB2 NM_024675.3 601355 Partner and localizer of BRCA2 4069 13
PTEN NM_000314.4 601728 Phosphatase and tensin 5572 9
RAD50 NM_005732.3 604040 DNA repair protein RAD50 homolog 6597 25
RAD51 NM_001164269.1 179617 DNA repair protein RAD51A homolog 2147 10
STK11 NM_000455.4 602216 Serine/threonine protein kinase 11 3286 10
TGFB1 NM_000660.4 190180 Transforming growth factor beta1 2217 7
TP53 NM_000546.5 191170 Tumor protein p53 2591 11