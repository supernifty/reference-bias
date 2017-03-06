#!/bin/bash
module load python-gcc/2.7.5
source ./venv/bin/activate
module load bwa-intel/0.7.12
module load samtools-intel/1.2
module load bedtools-intel/2.17.0
module load mauve/2.4.0
#sbatch ./xpipe-china-16.sh
#sbatch ./xpipe-korea-22.sh
#sbatch ./xpipe-hg19-22.sh
#sbatch ./xpipe-hg18-22.sh
#sbatch ./xpipe-china-22.sh

# hg19
#python extract_fasta.py chr1 < ./data/human/hg19.fa > ./data/human/hg19.1.fa
#python extract_fasta.py chr2 < ./data/human/hg19.fa > ./data/human/hg19.2.fa
#python extract_fasta.py chr3 < ./data/human/hg19.fa > ./data/human/hg19.3.fa
#python extract_fasta.py chr4 < ./data/human/hg19.fa > ./data/human/hg19.4.fa
#python extract_fasta.py chr5 < ./data/human/hg19.fa > ./data/human/hg19.5.fa
#python extract_fasta.py chr6 < ./data/human/hg19.fa > ./data/human/hg19.6.fa
#python extract_fasta.py chr7 < ./data/human/hg19.fa > ./data/human/hg19.7.fa
#python extract_fasta.py chr8 < ./data/human/hg19.fa > ./data/human/hg19.8.fa
#python extract_fasta.py chr9 < ./data/human/hg19.fa > ./data/human/hg19.9.fa
#python extract_fasta.py chr10 < ./data/human/hg19.fa > ./data/human/hg19.10.fa
#python extract_fasta.py chr11 < ./data/human/hg19.fa > ./data/human/hg19.11.fa
#python extract_fasta.py chr12 < ./data/human/hg19.fa > ./data/human/hg19.12.fa
#python extract_fasta.py chr13 < ./data/human/hg19.fa > ./data/human/hg19.13.fa
#python extract_fasta.py chr14 < ./data/human/hg19.fa > ./data/human/hg19.14.fa
#python extract_fasta.py chr15 < ./data/human/hg19.fa > ./data/human/hg19.15.fa
#python extract_fasta.py chr16 < ./data/human/hg19.fa > ./data/human/hg19.16.fa
#python extract_fasta.py chr17 < ./data/human/hg19.fa > ./data/human/hg19.17.fa
#python extract_fasta.py chr18 < ./data/human/hg19.fa > ./data/human/hg19.18.fa
#python extract_fasta.py chr19 < ./data/human/hg19.fa > ./data/human/hg19.19.fa
#python extract_fasta.py chr20 < ./data/human/hg19.fa > ./data/human/hg19.20.fa
#python extract_fasta.py chr21 < ./data/human/hg19.fa > ./data/human/hg19.21.fa
#python extract_fasta.py chr22 < ./data/human/hg19.fa > ./data/human/hg19.22.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-1.out ./data/human/hg19.1.fa ./data/human/hg38.1.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-2.out ./data/human/hg19.2.fa ./data/human/hg38.2.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-3.out ./data/human/hg19.3.fa ./data/human/hg38.3.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-4.out ./data/human/hg19.4.fa ./data/human/hg38.4.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-5.out ./data/human/hg19.5.fa ./data/human/hg38.5.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-6.out ./data/human/hg19.6.fa ./data/human/hg38.6.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-7.out ./data/human/hg19.7.fa ./data/human/hg38.7.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-8.out ./data/human/hg19.8.fa ./data/human/hg38.8.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-9.out ./data/human/hg19.9.fa ./data/human/hg38.9.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-10.out ./data/human/hg19.10.fa ./data/human/hg38.10.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-11.out ./data/human/hg19.11.fa ./data/human/hg38.11.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-12.out ./data/human/hg19.12.fa ./data/human/hg38.12.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-13.out ./data/human/hg19.13.fa ./data/human/hg38.13.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-14.out ./data/human/hg19.14.fa ./data/human/hg38.14.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-15.out ./data/human/hg19.15.fa ./data/human/hg38.15.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-16.out ./data/human/hg19.16.fa ./data/human/hg38.16.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-17.out ./data/human/hg19.17.fa ./data/human/hg38.17.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-18.out ./data/human/hg19.18.fa ./data/human/hg38.18.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-19.out ./data/human/hg19.19.fa ./data/human/hg38.19.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-20.out ./data/human/hg19.20.fa ./data/human/hg38.20.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-21.out ./data/human/hg19.21.fa ./data/human/hg38.21.fa
#progressiveMauve --output=./tmp_hg19/mauve-human-22.out ./data/human/hg19.22.fa ./data/human/hg38.22.fa

# hg18
#python extract_fasta.py chr1 < ./data/human/hg18.fa > ./data/human/hg18.1.fa
#python extract_fasta.py chr2 < ./data/human/hg18.fa > ./data/human/hg18.2.fa
#python extract_fasta.py chr3 < ./data/human/hg18.fa > ./data/human/hg18.3.fa
#python extract_fasta.py chr4 < ./data/human/hg18.fa > ./data/human/hg18.4.fa
#python extract_fasta.py chr5 < ./data/human/hg18.fa > ./data/human/hg18.5.fa
#python extract_fasta.py chr6 < ./data/human/hg18.fa > ./data/human/hg18.6.fa
#python extract_fasta.py chr7 < ./data/human/hg18.fa > ./data/human/hg18.7.fa
#python extract_fasta.py chr8 < ./data/human/hg18.fa > ./data/human/hg18.8.fa
#python extract_fasta.py chr9 < ./data/human/hg18.fa > ./data/human/hg18.9.fa
#python extract_fasta.py chr10 < ./data/human/hg18.fa > ./data/human/hg18.10.fa
#python extract_fasta.py chr11 < ./data/human/hg18.fa > ./data/human/hg18.11.fa
#python extract_fasta.py chr12 < ./data/human/hg18.fa > ./data/human/hg18.12.fa
#python extract_fasta.py chr13 < ./data/human/hg18.fa > ./data/human/hg18.13.fa
#python extract_fasta.py chr14 < ./data/human/hg18.fa > ./data/human/hg18.14.fa
#python extract_fasta.py chr15 < ./data/human/hg18.fa > ./data/human/hg18.15.fa
#python extract_fasta.py chr16 < ./data/human/hg18.fa > ./data/human/hg18.16.fa
#python extract_fasta.py chr17 < ./data/human/hg18.fa > ./data/human/hg18.17.fa
#python extract_fasta.py chr18 < ./data/human/hg18.fa > ./data/human/hg18.18.fa
#python extract_fasta.py chr19 < ./data/human/hg18.fa > ./data/human/hg18.19.fa
#python extract_fasta.py chr20 < ./data/human/hg18.fa > ./data/human/hg18.20.fa
#python extract_fasta.py chr21 < ./data/human/hg18.fa > ./data/human/hg18.21.fa
#python extract_fasta.py chr22 < ./data/human/hg18.fa > ./data/human/hg18.22.fa

#zcat ./data/human/china.fa | python extract_fasta.py chr16 > ./data/human/china-chr16.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr15 > ./data/human/china-chr15.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr14 > ./data/human/china-chr14.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr13 > ./data/human/china-chr13.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr12 > ./data/human/china-chr12.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr11 > ./data/human/china-chr11.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr10 > ./data/human/china-chr10.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr9 > ./data/human/china-chr9.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr8 > ./data/human/china-chr8.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr7 > ./data/human/china-chr7.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr6 > ./data/human/china-chr6.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr5 > ./data/human/china-chr5.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr4 > ./data/human/china-chr4.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr3 > ./data/human/china-chr3.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr2 > ./data/human/china-chr2.fa
#zcat ./data/human/china.fa | python extract_fasta.py chr1 > ./data/human/china-chr1.fa


#sbatch ./xpipe-template/xpipe-hg19.sh
sbatch ./xpipe-template/xpipe-hg18.sh
#sbatch ./xpipe-template/xpipe-china.sh
#sbatch ./xpipe-template/xpipe-china-self.sh

#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr1 > ./data/human/korea.1.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr2 > ./data/human/korea.2.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr3 > ./data/human/korea.3.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr4 > ./data/human/korea.4.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr5 > ./data/human/korea.5.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr6 > ./data/human/korea.6.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr7 > ./data/human/korea.7.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr8 > ./data/human/korea.8.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr9 > ./data/human/korea.9.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr10 > ./data/human/korea.10.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr11 > ./data/human/korea.11.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr12 > ./data/human/korea.12.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr13 > ./data/human/korea.13.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr14 > ./data/human/korea.14.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr15 > ./data/human/korea.15.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr16 > ./data/human/korea.16.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr17 > ./data/human/korea.17.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr18 > ./data/human/korea.18.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr19 > ./data/human/korea.19.fa
#zcat ./data/human/korea.fa.gz | python extract_fasta.py chr21 > ./data/human/korea.21.fa
#
#sbatch ./xpipe-template/xpipe-korea.sh
