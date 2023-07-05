#!/bin/bash

##### Data cleaning RNA-seq data ######

########################################
## Quality assessment 1
## project1/results/results_fastqc/fastqc_raw
########################################
## Assess read quality using FastQC

mkdir -p ~/project1/results/results_fastqc/fastqc_raw/

for filename in ~/project1/kadm-data/fastq/*fastq.gz ;
do
echo "Starting FastQC analysis of" $filename

fastqc $filename \
       -o ~/project1/results/results_fastqc/fastqc_raw/ \
       -t 15
done

########################################
## Adapter removal
## project1/data/fastq_trim
########################################
## Check if adapters exist

AdapterRemoval --identify-adapters \
    --file1 ~/project1/kadm-data/fastq/AM10_AM_Media_S1_R1_001.fastq.gz \
    --file2 ~/project1/kadm-data/fastq/AM10_AM_Media_S1_R2_001.fastq.gz \
    --threads 15

## Remove adapters 
## Remove reads with > 1 ambiguous base
## Trim ends until reach base with quality 30+
## Remove reads < 15 bp

mkdir -p ~/project1/data/fastq_trim

paste <(ls ~/project1/kadm-data/fastq/*R1_001.fastq.gz) \
      <(ls ~/project1/kadm-data/fastq/*R2_001.fastq.gz) |

while read file1 file2;
do
  name1=$(paste -d '\0' \
            <(echo 'project1/data/fastq_trim/') \
            <(awk -F'[_]S' '{print $1}' <(basename $file1)))
  
  AdapterRemoval --file1 $file1 --file2 $file2 \
    --basename $name1 --gzip \
    --trim5p 9 --maxns 1 --minlength 15 \
    --trimqualities --minquality 30 \
    --threads 15
done
      
########################################
## Quality assessment 2
## project1/results/results_fastqc/fastqc_trim
########################################
## Check if adapters still exist

AdapterRemoval --identify-adapters \
    --file1 ~/project1/data/fastq_trim/AM10_AM_Media.pair1.truncated.gz \
    --file2 ~/project1/data/fastq_trim/AM10_AM_Media.pair2.truncated.gz \
    --threads 15
    
## Assess trimmed read quality using FastQC.

mkdir -p ~/project1/results/results_fastqc/fastqc_trim/

for filename in ~/project1/data/fastq_trim/*pair[12].truncated.gz ;
do
echo "Starting FastQC analysis of" $filename
fastqc $filename \
       -o ~/project1/results/results_fastqc/fastqc_trim/ \
       -t 15
done
    
########################################
## Alignment
## project1/data/bam
########################################
## Get ref data files
mkdir -p ~/project1/data/STARindex
mkdir -p ~/project1/data/STARref
cd ~/project1/data/STARref

sudo curl -O ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz

sudo curl -O ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

cd

## Make genome index
STAR --runMode genomeGenerate \
     --genomeDir ~/project1/data/STARindex \
     --genomeFastaFiles ~/project1/data/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile ~/project1/data/STARref/Homo_sapiens.GRCh38.99.gtf \
     --sjdbOverhang 99 \
     --runThreadN 15
     
## Align with STAR
mkdir -p ~/project1/data/bam/

paste <(ls ~/project1/data/fastq_trim/*pair1.truncated.gz) \
      <(ls ~/project1/data/fastq_trim/*pair2.truncated.gz) |

while read file1 file2;
do
    echo "Aligning" $(basename  -- "$file1");
    
    name=$(paste -d '\0' \
            <(echo 'project1/data/bam/') \
            <(awk -F'[.]pair' '{print $1}' <(basename $file1)) \
            <(echo '_'))
    
    STAR --genomeDir project1/data/STARindex \
         --readFilesIn $file1 $file2 \
         --readFilesCommand zcat \
         --outFileNamePrefix $name \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 15 \
         --runRNGseed 8756
done

########################################
## Assess aligned reads
## project1/results/results_cleaning
########################################
## Get refFlat genome for Picard
mkdir -p ~/project1/data/PICARDref
cd ~/project1/data/PICARDref
sudo curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
## Remove chr in chromosome name to match ensembl alignment
sed 's/chr//' refFlat.txt > refFlat.ensembl.txt
cd

## median CV of gene model coverage

mkdir -p ~/project1/results/results_cleaning/

for bam_file in ~/project1/data/bam/*sortedByCoord.out.bam ;
do
    java -XX:ParallelGCThreads=15 \
        -jar ~/project1/applications/anaconda3/share/picard-2.22.0-0/picard.jar \
        CollectRnaSeqMetrics \
        REF_FLAT=~/project1/data/PICARDref/refFlat.ensembl.txt \
        INPUT=$bam_file  \
        OUTPUT=~/project1/results/results_cleaning/temp.tsv \
        ASSUME_SORTED=true \
        STRAND_SPECIFICITY=NONE \
        MINIMUM_LENGTH=500
    
    #Append results
    echo $bam_file >> ~/project1/results/results_cleaning/bam.metrics.tsv
    cat ~/project1/results/results_cleaning/temp.tsv >> ~/project1/results/results_cleaning/bam.metrics.tsv
    #Remove this iteration
    rm ~/project1/results/results_cleaning/temp.tsv
done
        
## mapped_reads_w_dups aka alignment percentage

for bam_file in ~/project1/data/bam/*sortedByCoord.out.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> ~/project1/results/results_cleaning/summary.alignment.tsv
    
    samtools flagstat -@ 15 $bam_file \
    >> ~/project1/results/results_cleaning/summary.alignment.tsv
done

########################################
## Quality filter BAM
## project1/data/bam_filter
########################################
## -h: keep header
## -f 3: keeps paired reads where both mapped
## -F 1284: removes unmapped reads, non-primary alignments, and PCR duplicates
## -q 30: min MAPQ of 30

## Paired reads
mkdir -p project1/data/bam_filter_paired/

for bam_file in project1/data/bam/*sortedByCoord.out.bam ;
do
  filename1=$(paste -d '\0' \
            <(echo 'project1/data/bam_filter_paired/') \
            <(awk -F'[_]Aligned' '{print $1}' <(basename $bam_file)) \
            <(echo '_filter_paired.bam'))
  
  echo "Filtering" $bam_file          
  samtools view $bam_file \
      -h -f 3 -F 1284 -q 30 \
      -@ 15 \
      > $filename1
done

## alignment percentage

for bam_file in project1/data/bam_filter_paired/*filter_paired.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> project1/results/results_cleaning/summary.align.filter.paired.tsv
    
    samtools flagstat -@ 15 $bam_file \
    >> project1/results/results_cleaning/summary.align.filter.paired.tsv
done

########################################
## Count reads in genes
## project1/data/counts
########################################
## List all possible features to count in annotation file
#cut -f 3 project1/data/STARref/Homo_sapiens.GRCh38.99.gtf | sort | uniq
mkdir -p project1/data/counts

## Count reads in genes
featureCounts -T 14 -g gene_id -t exon -p \
  -a project1/data/STARref/Homo_sapiens.GRCh38.99.gtf \
  -o project1/data/counts/AM.MDM.featurecounts.paired.tsv \
  project1/data/bam_filter_paired/*filter_paired.bam

################# END ##################