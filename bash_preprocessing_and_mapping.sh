#!/bin/bash

mkdir demult

sbatch --time 20:00:00 --mem 32G --wrap "/Users/simeon01/applications/demultiplexer.rhel/demuxFQ -c -d -i -e -t 1 -r 0.01 -R -l 9 -o ./demult -b ./demult/SLX-29000.lostreads.r_1.fq.gz -s SLX-29000_summary_demuxFQ_r1.txt SLX-29000_demuxFQ_r1.txt SLX-29000.UnspecifiedIndex.HXXXXXXX.s_1.r_1.fq.gz"

sbatch --mem 32G --time 20:00:00 --wrap "/Users/simeon01/applications/demultiplexer.rhel/demuxFQ -c -d -i -e -t 1 -r 0.01 -R -l 9 -o ./demult -b ./demult/SLX-29000.lostreads.r_2.fq.gz -s SLX-29000_summary_demuxFQ_r2.txt SLX-29000_demuxFQ_r2.txt SLX-29000.UnspecifiedIndex.HXXXXXXX.s_1.r_2.fq.gz"





#!/bin/bash

for file in *.fq.gz
do
sbatch --time 02:00:00 -o %j.out -e %j.err --mem 10G --wrap "fastqc $file"
done





#!/bin/bash

mkdir trimmed

for file in *R1.fq.gz
do
fq1=$file
fq2=${fq1/R1/R2}

sbatch --time 03:00:00 -o %j.out -e %j.err --mem 30G --wrap "cutadapt -q 20  -o trimmed/${fq1%%.fq.gz}.trimmed.fq.gz -p trimmed/${fq2%%.fq.gz}.trimmed.fq.gz $fq1 $fq2"
done




#!/bin/bash

mkdir -p out.$1

#Split fastq files
seqkit split2 -1 $1_R1.fq.gz -2 $1_R2.fq.gz -s 4000000 -O out.$1 -f -e .gz

cd out.$1

#Reference genome
REF='/Users/dhir01/References/Genomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna'

#Read alignment (bwa mem)
for fq1 in *R1.part*.fq.gz; do fq2=${fq1/R1/R2}; bname=${fq1%.fq.gz}; echo $fq1,$fq2,$bname; sbatch --time 10:00:00 --mem 90G -e %j.err --wrap "bwa mem -SP5M -t 18 $REF $fq1 $fq2 >$bname.bam";done





#!/bin/bash

out=$1
mapq=$2
header=$1'_R1.part_001.bam'


#merge split files
samtools merge -o $out.bam -h $header *_R1.part*bam


samtools view -h $out.bam | pairtools parse -c /Users/dhir01/References/Genomes/GRCh38/hg38.chr.genome --min-mapq $mapq --assembly hg38 -o $out.parsed.pairsam
pairtools sort -o $out.sorted.pairsam $out.parsed.pairsam
pairtools dedup --mark-dups -o $out.deduped.pairsam $out.sorted.pairsam --output-stats $out.stats

pairtools select '(pair_type == "UU" or pair_type == "UR" or pair_type == "RU")' --output $out.filtered.pairsam $out.deduped.pairsam
pairtools split --output-pairs $out.nodups.pairs.gz $out.filtered.pairsam

zcat $out.nodups.pairs.gz | grep -v '#'  | perl -F'\t' -lane 'print join "\t", @F[0,1,2,5,3,4,6]' >$out.Valid.pairs
gzip $out.Valid.pairs  
