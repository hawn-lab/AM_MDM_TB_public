#!/bin/bash

##### Data cleaning RNA-seq data ######
mkdir ~/project/data2
sudo chmod 777 -R ~/project/data2
s3fs kadm-results-am ~/project/data2 -o passwd_file=~/.passwd-s3fs \
    -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007
cd

########################################
## Alignment
########################################
## Get ref data files
mkdir -p ~/project/Mtb_index
cd ~/project/Mtb_index

tar -C . -xvf ~/project/ref/Mtb.H37Rv.NC_000962.3/genome_assemblies_genome_gtf.tar
gunzip ncbi-genomes-2021-03-01/GCF_000195955.2_ASM19595v2_genomic.gtf.gz
    
tar -C . -xvf ~/project/ref/Mtb.H37Rv.NC_000962.3/genome_assemblies_genome_fasta.tar
gunzip ncbi-genomes-2021-03-01/GCF_000195955.2_ASM19595v2_genomic.fna.gz

aws s3 sync ncbi-genomes-2021-03-01/ s3://kadm-ref/Mtb.H37Rv.NC_000962.3/

## Make genome index
STAR --runMode genomeGenerate \
     --genomeDir ~/project/Mtb_index \
     --genomeFastaFiles ~/project/ref/Mtb.H37Rv.NC_000962.3/GCF_000195955.2_ASM19595v2_genomic.fna \
     --sjdbGTFfile ~/project/ref/Mtb.H37Rv.NC_000962.3/GCF_000195955.2_ASM19595v2_genomic.gtf \
     --sjdbOverhang 99 \
     --runThreadN 60

cp ./Log.out ~/project/Mtb_index/
aws s3 sync ~/project/Mtb_index s3://kadm-ref/Mtb.H37Rv.NC_000962.3/

## Align with STAR
mkdir -p ~/project/results.mtb/bam/

#regex is hard. Run sample 69
STAR --genomeDir project/ref/Mtb.H37Rv.NC_000962.3/ \
         --readFilesIn ~/project/data2/fastq_trim/69_AM18_Cryo_MTB.pair1.truncated.gz \
         ~/project/data2/fastq_trim/69_AM18_Cryo_MTB.pair2.truncated.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix project/results.mtb/bam/69_AM18_Cryo_MTB_ \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 20 \
         --runRNGseed 8756
#then run samples 70-74
paste <(ls ~/project/data2/fastq_trim/7[0-4]*pair1.truncated.gz) \
      <(ls ~/project/data2/fastq_trim/7[0-4]*pair2.truncated.gz) |

while read file1 file2;
do
    echo "Aligning" $(basename  -- "$file1");
    
    name=$(paste -d '\0' \
            <(echo 'project/results.mtb/bam/') \
            <(awk -F'[.]pair' '{print $1}' <(basename $file1)) \
            <(echo '_'))
    
    STAR --genomeDir project/ref/Mtb.H37Rv.NC_000962.3/ \
         --readFilesIn $file1 $file2 \
         --readFilesCommand zcat \
         --outFileNamePrefix $name \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 20 \
         --runRNGseed 8756
done

########################################
## Quality filter BAM
########################################

## -h: keep header
## -f 3: keeps paired reads where both mapped
## -F 1284: removes unmapped reads, non-primary alignments, and PCR duplicates
## -q 30: min MAPQ of 30

## Paired reads
mkdir -p project/results.mtb/bam_filter_paired/

for bam_file in project/results.mtb/bam/*sortedByCoord.out.bam ;
do
  filename1=$(paste -d '\0' \
            <(echo 'project/results.mtb/bam_filter_paired/') \
            <(awk -F'[_]Aligned' '{print $1}' <(basename $bam_file)) \
            <(echo '_filter_paired.bam'))
  
  echo "Filtering" $bam_file          
  samtools view $bam_file \
      -h -f 3 -F 1284 -q 30 \
      -@ 40 \
      > $filename1
done

########################################
## Assess aligned reads
########################################
## median CV of gene model coverage

mkdir -p ~/project/results.mtb/results_cleaning

for bam_file in ~/project/results.mtb/bam/*sortedByCoord.out.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> ~/project/results.mtb/results_cleaning/summary.alignment.tsv
    
    samtools flagstat -@ 60 $bam_file \
    >> ~/project/results.mtb/results_cleaning/summary.alignment.tsv
done

###

for bam_file in ~/project/results.mtb/bam_filter_paired/*filter_paired.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> ~/project/results.mtb/results_cleaning/summary.alignment.tsv
    
    samtools flagstat -@ 60 $bam_file \
    >> ~/project/results.mtb/results_cleaning/summary.alignment.tsv
done

