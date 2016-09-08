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

