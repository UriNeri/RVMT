#!/bin/bash
#hostname

THREADS=10
work_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/Valerian/IDFT_blastn_ICTV/"
input="/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/Contigs/IDFT.fna"
ini_name=$(basename $input ".fna")

cd $work_dir

pid=10
makeblastdb -in RefSeqs2Get.fas -dbtype 'nucl' -out ICTVdb/ICTVdb

blastn -task dc-megablast -evalue 0.001 -perc_identity 66 -word_size 12 -max_target_seqs 10000000 -num_threads $THREADS -db ICTVdb/ICTVdb -out /urigo/urineri/DCblastn_IDFT_vs_V302filt.tsv -outfmt "6 qseqid sseqid pident sstart send qstart qend slen qlen length evalue bitscore"
blastn -task dc-megablast -culling_limit 10 -perc_identity $pid -query "$input" -db ICTVdb/ICTVdb -out blastn_IDFT_vs_ICTV_out.tsv -num_threads $THREADS -outfmt "6 qseqid sseqid qcovs nident positive pident qlen slen qstart qend sstart send length evalue bitscore gaps"
