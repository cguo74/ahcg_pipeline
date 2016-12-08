# Fall2016 BIOL-8803F
# *Project1 Running ahcg_pipeline*

###Variant calling pipeline for genomic data analysis###

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
```

## Copy the pipeline
```{sh}
search ahcg_pipeline in github then folk it
https://github.com/shashidhar22/ahcg_pipeline
```

## Changing a remote's URL
```{sh}
instruction: https://help.github.com/articles/changing-a-remote-s-url/

git remote set-url https://github.com/cguo74/ahcg_pipeline.git
git remote -v

below means that the Changing a remote's URL wasn’t success
origin	https://github.com/shashidhar22/ahcg_pipeline (fetch)
origin	https://github.com/shashidhar22/ahcg_pipeline (push)

we add origin
git remote set-url origin https://github.com/cguo74/ahcg_pipeline.git
git remote -v

below means that the Changing a remote's URL was success
origin	https://github.com/cguo74/ahcg_pipeline.git (fetch)
origin	https://github.com/cguo74/ahcg_pipeline.git (push)
```
## Download emacs(easy to use than vim)
```{sh}
sudo apt-get update
sudo apt-get install emacs24
```

## Ignore the unnecessary folder, add any folder under current position
```{sh}
basespace@ubuntu:~/ahcg_pipeline$ ls -a 
emacs .gitignore
```


## Fasta index using Samtools faidx

###Download samtools will produce hg19.fa.fai###
```{sh}
sudo apt-get install samtools
instructions: open the page: http://www.htslib.org/doc/samtools.html
```
####final commend####
```{sh}
samtools faidx hg19.fa
```


## Genome dict file using picard will produce hg19.dict
###download java and picard###
```{sh}
sudo apt-get install openjdk-8-jre-headless
sudo apt-get install picard
```

###Also download picard jar here https://github.com/broadinstitute/picard/releases/tag/2.6.0###
```{sh}
java -jar picard.jar
instructions: https://broadinstitute.github.io/picard/command-line overview.html#CreateSequenceDictionary
```
####final commend####
```{sh}
java -jar picard.jar CreateSequenceDictionary R=hg19.fa O=reference.dict
```

###Something need to correct before running.In the folder named wantignore I have all the files:patient genome-test_r1.fastq; references genome after bowtie processed- hg19.1.bt2; 
###In resources/genome, I have hg19.fa which is original file, hg19.fa.fai which generated from samtool, hg19.dict generated from picard###
```{sh}
mv reference.dict ../resources/genome
mv hg19.fa.fai ../resources/genome
gunzip ./resources/dbsnp/dbsnp_138.hg19.vcf.gz
```
## Try run the pipeline
###To access help use the following command:###
```{sh}
python3 ahcg_pipeline.py -h
```

##Final commend for first project
```{sh}
python3 ahcg_pipeline.py -t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b ./lib/bowtie2-2.2.9/bowtie2 -p ./lib/picard.jar -g ./lib/GenomeAnalysisTK.jar -i ./wantignore/test_r1.fastq ./wantignore/test_r2.fastq -w ./wantignore/hg19 -d ./resources/dbsnp/dbsnp_138.hg19.vcf -r ./resources/genome/hg19.fa -a ./lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o ./hw1
```
#Update an existing file to a GitHub repository
###Copy the file that needed to be updated/uploaded to GitHub into the local directory that was created when you cloned the repository###
###Change the current working directory to local repository###

###Stage the file for commit to your local repository using the following command###
```{sh}
$ git add .
This command adds / stages all of the files in the current directory. This is for convenience, and can still be used if you have certain files you don't want to add by using a .gitignore
```
###Commit the file that you've staged in local repository###
```{sh}
$ git commit -m "Notes related with changes"

This command stores the current contents of the index in a new commit along with a log message from the user describing the changes
```
###Push the chagnes in local repository to Github remote repository###
```{sh}
$ git push remote_name branch_name
For example, use $ git push origin master
```

# *Project2 Extracting reads mapping to BRCA1 from NA12878 HiSeq Exome dataset*

## Download the reference genome file
```{sh}
$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
```

## Extract the genome coordinates for BRCA1
```{sh}
$ grep -w "NM_007294" hg19_refGene.txt > NM_007294.txt
  Note: NM_007294 is the popular transcript of BRCA1
```

## Extract chrom, chromStart, chromEnd information from NM_007294 and convert it to bed file format
```{sh}
$ awk '{print $3, $5-1000, $6+1000}' NM_007294.txt > NM_007294_gene_coordinates.bed
```

## Download the NA12878 HiSeq Exome dataset BAM files from the link

## Using samtools to subset the bam file to regions corresponding to BRCA1
```{sh}
$ samtools view -L <bed file> -b -o < outout bam file > < input bam file >
    Note: -b just specifies that the output needs to be a bam file.
    Example code: $ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_1.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam
$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_2.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam
$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_3.bam project.NIST_NIST7086_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam
$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_4.bam project.NIST_NIST7086_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam
```
## Merge bam files
```{sh}
$ samtools merge finalBamFile.bam *.bam
```

## Using bedtools to convert the bam file to a fastq file
```{sh}
$ samtools sort -n aln.bam aln.qsort.bam
$ bedtools bamtofastq -i <bam file> -fq < fastq r1> -fq2 < fastq r2>
Note:From the BRCA1 bam file we now extract the reads aligning to the region using bedtools
Note:If your BAM alignments are from paired-end sequence data, one can use the -fq2 option to create two distinct FASTQ output files, one for end 1 and one for end 2. When using this option, it is required that the BAM file is sorted/grouped by the read name. This keeps the resulting records in the two output FASTQ files in the same order.
```
## Picard genome dictionary
```{sh}
$ java -jar picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict   Note: The reference file could be download from link
```
## Samtools fasta index
```{sh}
$ samtools faidx <path to ref.fa>
```

## Bowtie index
```{sh}
$ bowtie2-build <path to ref.fa> <output prefix>
```

## Run the ahcg_pipeline.py to get BRCA1 variants
```{sh}
$ python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i lab/BRCA1_bam/merge_BRCA1_end1.fq lab/BRCA1_bam/merge_BRCA1_end2.fq -w lab/chr17 -d dbsnp/dbsnp_138.hg19.vcf -r lab/chr17.fa -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o output_merge_BRC1
```

# *Project3 Compare ahcg_pipeline_vcf with golden standard vcf file*

## Download the BED file, which contains the Nextera Rapid Capture Expanded Exome targeted regions' genomic coordinates
```{sh}
$ wget http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_expandedexome_targetedregions.bed
```
## Download the vcf file running from ahcg_pipeline
```{sh}
$ wget http://vannberg.biology.gatech.edu/data/variants.vcf
```

## Download the golden standard vcf file
```{sh}
$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/10XGenomics_calls_08142015/NA12878_phased_variants.vcf.gz
```
### Use bedtools to cut the Nextera Rapid Capture Expanded Exome targeted regions from golden standard vcf file and ahcg_pipeline vcf file
```{sh}
$ bedtools intersect -wa -a <vcf file> -b <bed file>
```

### Use bedtools to find intersect calls between ahcg_pipeline vcf and golden standard vcf
```{sh}
$ bedtools intersect -header -a <vcf_golden_standard> -b <vcf_pipeline>
```

### Use bedtools subtract command to find calls that are unique for gold_standard vcf
```{sh}
$ bedtools substract -a <vcf_golden_standard> -b <vcf_intersect>
```

### Use bedtools subtract commandto find calls that are unique for ahcg_pipeline vcf
```{sh}
$ bedtools substract -a <vcf_pipeline> -b <vcf_intersect>
```

### Use bedtools jaccard to find the number of TP, FP, FN variant calls
```{sh}
$ bedtools jaccard -a <vcf_intersect> -b <vcf_pipeline>
$ bedtools jaccard -a <vcf_intersect> -b <vcf_golden_standard>
```

# *Project4 Extract chromosome coordinates given NCBI_Accession number*

### Download the reference genome file
```{sh}
$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
```

### Download the gene list of breast cancer biomarkers. link
```{sh}
    Note: Extract the first two columns and exported as a txt file.
    Note: The first column is the gene symbols, the second column is the NCBI reference.
```

### Run the python script to get the genome coordinates for the NCBI_Accession list
```{sh}
$ python get_coordinates_from_gene.py
```

# *Project5 Compare ahcg_pipeline_vcf with GIAB_vcf for BRC_OC_genes*

### Download the vcf file running from ahcg_pipeline

### Download the golden vcf file running from GIAB.
```{sh}
$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz
```

### Make the BED file for the gene related with breast cancer(BRC) and ovarian cancer(OC)
```{sh}
1. Gene list from color genomics
2. Gene list from otogenetics
3. Create the gene list with gene symbols and associated transcript accession numbers.
    Note: The first column is the gene symbols, the second column is its associated transcript's NCBI reference number.
4. Download the reference genome file
$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
5. Run the python script to get the genome coordinates for the gene list
$ python get_coordinates_from_gene.py
6. Add 20bp flank region to both the chrom_star and chrom_end of each gene's chromosome coordinates.
$ python 20bp_flank_both_ends.py
```
### Use bedtools to extract variants from both ahcg_pipeline vcf files and GIAB golden vcf, given the bed file of breast cancer and ovarian cancer related genes' chromosome coordinates
```{sh}
$ bedtools intersect -wa -header -a <vcf file> -b <bed file>
```

### Use bedtools to find intersect variants between vcf(from ahcg_pipeline) and golden standard vcf
```{sh}
$ bedtools intersect -a <vcf_golden_standard> -b <vcf_ahcg_pipeline>
```

# *Project6 Apply GATK-variantRecalibrator on raw vcf file*

### Download the vcf file running from ahcg_pipeline

### Download the required known variant calls for VQSR to build model
```{sh}
  Note: Creating a tabix indexed vcf file
$ gunzip file.gz
$ bgzip -c file.vcf > file.vcf.gz
$ tabix -p vcf file.vcf.gz
```
### VariantRecalibrator - build a recalibration model to score variant quality for filtering purposes
```{sh}
$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar \ -T VariantRecalibrator \ -R ref_genome/hg19.fa \ -input vcf_files_not_upload/NA12878_variants_ahcg.vcf \ -resource:hapmap,known=false,training=true,truth=true,prior=15.0 variantRecali/hapmap_3.3.hg19.sites.vcf \ -resource:omni,known=false,training=true,truth=false,prior=12.0 variantRecali/1000G_omni2.5.hg19.sites.vcf \ -resource:1000G,known=false,training=true,truth=false,prior=10.0 variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf \ -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp/dbsnp_138.hg19.vcf -an DP \ -an QD \ -an FS \ -an SOR \ -an MQRankSum \ -an ReadPosRankSum \ -mode SNP \ -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \ -recalFile recalibrate_SNP.recal \ -tranchesFile recalibrate_SNP.tranches
```
### ApplyRecalibration - Apply a score cutoff to filter variants based on a recalibration table
```{sh}
$ java -jar lib/GenomeAnalysisTK.jar \ -T ApplyRecalibration \ -R ref_genome/hg19.fa \ -input vcf_files_not_upload/NA12878_variants_ahcg.vcf \ -mode SNP \ --ts_filter_level 99.0 \ -recalFile recalibrate_SNP.recal \ -tranchesFile recalibrate_SNP.tranches \ -o recalibrated_snps_raw_indels.vcf
```
# *First Presentation*
## Calculate read depth based on alignment file(BRCA1)

### Extract BRCA1 gene chromosome coordinates from "brc_oc_gene_list_bed_add_20.txt"
```{sh}
$ grep 'NM_007294' brc_oc_gene_list_bed_add_20.txt > brca1.bed
```
### Extract brca1 alignments
```{sh}
$ samtools view -L brca1.bed project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam -b > na12878.brca1.bam
Note: -L: only output alignments overlapping in the input bed file
Note: -b: output alignments in the bam format
```

### Computes and summarize the depth for brca1
```{sh}
$ bedtools genomecov -ibam na12878.brca1.bam -bga > na12878.brca1.bga.bed
Note: -ibam: BAM file as input for coverage.
Note: -bga: Reporting genome coverage for all positions in BEDGRAPH format.
Note: This is BedGraph format: chrom chromStart chromEnd dataValue
```
### Intersection between two bed files
```{sh}
$ bedtools intersect -split -a na12878.brca1.bga.bed -b brca1.bed -bed > brca1.final.bed
Note: when using the -split option, only the exon overlaps are reported
```
### Use the read_depth_calculation.py to get the read depth result
```{sh}
$ python read_depth_calculation.py
```

## Annotate the vcf file with pathogenic information(BRCA1)

### Download the BRCA variant pathogenic annotation information file
```{sh}
$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA1_brca_exchange_variants.csv
$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA2_brca_exchange_variants.csv
```
### Use the pathgenic_annotation_add_dp.py script to add pathogenic information to the variants
```{sh}
$ python pathgenic_annotation_add_dp.py
```

## Automate the pipeline to add read depth information and pathogenic information for the variants
```{sh}
Run the ahcg_pipeline with fastq file and create the vcf file.
Run GATK genome recalibrator on the result vcf file.
Use get_coordinates_from_gene.py script to create the bed file for the chromosome coordinates of given genes of interests.
Use bedtools to get the variants for the given bed file from the recalibrated vcf file.
Use the pathgenic_annotation_add_dp.py script to add pathogenic information to the variants.
Use samtools to extract alignments given bed file and then use bedtools to compute and summarize the depth given bed file.
Use read_depth_calculation.py script to calculate the reads depth based on bed file.
Use the report_depth_pathogenic.py to integrate the reads depth information with the pathogenic information for the variants.
```

# *Final Project*

## Download bam file of patient1
```{sh}
wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai
```
## Variant calling analysis using GATK
```{sh}
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T HaplotypeCaller -R input_files/ref_genome/hg19.fa -I Patient1_RG_MD_IR_BQ.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o p1_raw_variants.vcf
```{sh}

## Perform variants recalibration analysis using GATK-VariantRecalibrator to get more accurate results
```{sh}
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R input_files
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R input_files/ref_genome/hg19.fa -input p1_raw_variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile p1_recalibrate_SNP.recal -tranchesFile p1_recalibrate_SNP.tranches -o p1_recalibrated_snps_raw_indels.vcf
```
## Create dcm-related 6 genes' genome position bed file
 ```{sh}
 scripts/grep_gen_list.pl
 ```

## Extract the dcm-related 6 genes' genomic variants fro, getting intersect between bed file and vcf
```{sh}
bedtools intersect -wa -header -a p1_recalibrated_snps_raw_indels.vcf -b scripts/temp/exom_list.bed > patient1_dcm_final.vcf
```

## Calculat the reads depth information for DCM genes

### Obtain depth for each gene
 ```{sh}
 samtools view -L scripts/temp/exom_list.bed Patient1_RG_MD_IR_BQ.bam -b > patient1.dcm.bam
 ```
### Reporting genome coverage
 ```{sh}
bedtools genomecov -ibam patient1.dcm.bam -bga > patient1.dcm.bed
 ```
### Read coverage for each gene
```{sh}
bedtools intersect -split -a patient1.dcm.bed -b scripts/temp/exom_list.bed -bed > scripts/temp/patient1.final.bed
python scripts/depth_for_each_site.py
scripts/label.pl scripts/temp/patient1.final.full.txt patient1.final.full.labeled.txt
```
#### Python script does:
##### 1.Took intersect between vcf and bed files to get variants

##### 2.Matched the variants with depth information and exon position from final txt file

##### 3. Use the depths and exon position compare to cutoff, and calculate (cutoff - depth)

##### 4. Save the result to a vcf last column named as "P1_afterinter.vcfnew.vcf" ect


## Coverage plotting
```{sh}
cp patient1.final.full.labeled.txt scripts/txt/p1.final.txt
R CMD BATCH '--args scripts/txt/p1.final.txt' scripts/plots.R
R version 3.14 with package ggplot2 needed
```
## Vcf analysis reporting
```{sh}
cp patient1_dcm_final.vcf scripts/vcf/P1_afterinter.vcf
python scripts/min.py
cp scripts/vcf/*new.vcf SNP_report/
```


