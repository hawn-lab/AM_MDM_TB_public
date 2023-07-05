#!/bin/bash

##### Data cleaning RNA-seq data ######

s3fs kadm-data.am ~/project/data.am -o passwd_file=~/.passwd-s3fs \
    -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007
cd

########################################
## Quality assessment 1
## project/results.am/results_fastqc/fastqc_raw
########################################
## Assess read quality using FastQC

mkdir -p ~/project/results.am/results_fastqc/fastqc_raw/

for filename in ~/project/data.am/*fastq.gz ;
do
echo "Starting FastQC analysis of" $filename

fastqc $filename \
       -o ~/project/results.am/results_fastqc/fastqc_raw/ \
       -t 90
done

aws s3 sync ~/project/results.am/ s3://kadm-results-am

########################################
## Adapter removal
## project/data/fastq_trim
########################################
## Remove adapters 
## Remove reads with > 1 ambiguous base
## Trim ends until reach base with quality 30+
## Remove reads < 15 bp

mkdir -p ~/project/results.am/fastq_trim

paste <(ls ~/project/data.am/*R1_001.fastq.gz) \
      <(ls ~/project/data.am/*R2_001.fastq.gz) |

while read file1 file2;
do
  name1=$(paste -d '\0' \
            <(echo 'project/results.am/fastq_trim/') \
            <(awk -F'[_]S' '{print $1}' <(basename $file1)))
  
  AdapterRemoval --file1 $file1 --file2 $file2 \
    --basename $name1 --gzip \
    --trim5p 7 --maxns 1 --minlength 15 \
    --trimqualities --minquality 30 \
    --threads 80
done
      
########################################
## Quality assessment 2
## project/results.am/results_fastqc/fastqc_trim
########################################
## Assess trimmed read quality using FastQC.

mkdir -p ~/project/results.am/results_fastqc/fastqc_trim/

for filename in ~/project/results.am/fastq_trim/*pair[12].truncated.gz ;
do
echo "Starting FastQC analysis of" $filename
fastqc $filename \
       -o ~/project/results.am/results_fastqc/fastqc_trim/ \
       -t 80
done
    
########################################
## Alignment
## project/data/bam
########################################
## Get ref data files
mkdir -p ~/project/results.am/STARindex
mkdir -p ~/project/results.am/STARref
cd ~/project/results.am/STARref

sudo curl -O ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz
gunzip Homo_sapiens.GRCh38.102.gtf.gz

sudo curl -O ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

cd

## Make genome index
STAR --runMode genomeGenerate \
     --genomeDir ~/project/results.am/STARindex \
     --genomeFastaFiles ~/project/results.am/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile ~/project/results.am/STARref/Homo_sapiens.GRCh38.102.gtf \
     --sjdbOverhang 99 \
     --runThreadN 60

cp ./Log.out /home/ec2-user/project/results.am/STARindex/

## Align with STAR
mkdir -p ~/project/results.am/bam/

paste <(ls ~/project/results.am/fastq_trim/*pair1.truncated.gz) \
      <(ls ~/project/results.am/fastq_trim/*pair2.truncated.gz) |

while read file1 file2;
do
    echo "Aligning" $(basename  -- "$file1");
    
    name=$(paste -d '\0' \
            <(echo 'project/results.am/bam/') \
            <(awk -F'[.]pair' '{print $1}' <(basename $file1)) \
            <(echo '_'))
    
    STAR --genomeDir project/results.am/STARindex \
         --readFilesIn $file1 $file2 \
         --readFilesCommand zcat \
         --outFileNamePrefix $name \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 20 \
         --runRNGseed 8756
done

########################################
## Assess aligned reads
## project/results.am/results_cleaning
########################################
## median CV of gene model coverage

mkdir -p ~/project/results.am/results_cleaning/

for bam_file in ~/project/results.am/bam/*sortedByCoord.out.bam ;
do
    java -XX:ParallelGCThreads=50 \
        -jar project/apps/anaconda3/share/picard-2.24.2-0/picard.jar \
        CollectRnaSeqMetrics \
        REF_FLAT=~/project/ref/PICARDref/refFlat.ensembl.txt \
        I=$bam_file  \
        O=project/results.am/results_cleaning/temp.tsv \
        ASSUME_SORTED=true \
        STRAND_SPECIFICITY=NONE \
        MINIMUM_LENGTH=500

    #Append results
    echo $bam_file >> ~/project/results.am/results_cleaning/bam.metrics.tsv
    cat ~/project/results.am/results_cleaning/temp.tsv >> ~/project/results.am/results_cleaning/bam.metrics.tsv
    #Remove this iteration
    rm ~/project/results.am/results_cleaning/temp.tsv
done

java -XX:ParallelGCThreads=50 \
        -jar ~/project/apps/anaconda3/share/picard-2.24.2-0/picard.jar \
        CollectRnaSeqMetrics \
        REF_FLAT=~/project/ref/PICARDref/refFlat.ensembl.txt \
        I=~/project/data.am/bam/80_AM21_Cryo_IFNb_Aligned.sortedByCoord.out.bam  \
        O=~/project/results.am/results_cleaning/temp.tsv \
        ASSUME_SORTED=true \
        STRAND_SPECIFICITY=NONE \
        MINIMUM_LENGTH=500
CollectRnaSeqMetrics \
    -REF_FLAT ~/project/ref/PICARDref/refFlat.ensembl.txt
    -I ~/project/data.am/bam/80_AM21_Cryo_IFNb_Aligned.sortedByCoord.out.bam
    -O ~/project/results.am/results_cleaning/temp.tsv
    -ASSUME_SORTED true -STRAND_SPECIFICITY NONE -MINIMUM_LENGTH 500

## mapped_reads_w_dups aka alignment percentage

for bam_file in ~/project/results.am/bam/*sortedByCoord.out.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> ~/project/results.am/results_cleaning/summary.alignment.tsv
    
    samtools flagstat -@ 60 $bam_file \
    >> ~/project/results.am/results_cleaning/summary.alignment.tsv
done

########################################
## Quality filter BAM
## project/data/bam_filter
########################################
## -h: keep header
## -f 3: keeps paired reads where both mapped
## -F 1284: removes unmapped reads, non-primary alignments, and PCR duplicates
## -q 30: min MAPQ of 30

## Paired reads
mkdir -p project/results.am/bam_filter_paired/

for bam_file in project/results.am/bam/*sortedByCoord.out.bam ;
do
  filename1=$(paste -d '\0' \
            <(echo 'project/results.am/bam_filter_paired/') \
            <(awk -F'[_]Aligned' '{print $1}' <(basename $bam_file)) \
            <(echo '_filter_paired.bam'))
  
  echo "Filtering" $bam_file          
  samtools view $bam_file \
      -h -f 3 -F 1284 -q 30 \
      -@ 40 \
      > $filename1
done

## alignment percentage

for bam_file in project/results.am/bam_filter_paired/*filter_paired.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> project/results.am/results_cleaning/summary.align.filter.paired.tsv
    
    samtools flagstat -@ 60 $bam_file \
    >> project/results.am/results_cleaning/summary.align.filter.paired.tsv
done

########################################
## Count reads in genes
## project/data/counts
########################################
## List all possible features to count in annotation file
#cut -f 3 project/data/STARref/Homo_sapiens.GRCh38.99.gtf | sort | uniq
mkdir -p project/results.am/counts

## Count reads in genes
featureCounts -T 50 -g gene_id -t exon -p \
  -a project/ref/release102/STARref/Homo_sapiens.GRCh38.102.gtf \
  -o project/results.am/counts/AM.MDM.IFN.featurecounts.paired.tsv \
  project/results.am/bam_filter_paired/*filter_paired.bam

################# END ##################