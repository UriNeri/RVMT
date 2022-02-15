#!/bin/bash
#hostname

round=2.03
Attempt_id=3
iter=3
THREADS=7
DBpath="/media/uri/HDD1/uri/RNA_Vir_MTs/V2/TNseqs/AntiFam_6.0/AntiFam.hmm"
input_seqs="./uniq_Seqs4profiles_trimmed_Iter_3.faa"
input_name="Seqs4profiles_trimmed_Iter_3"
DB_name=AntiFam
cd ./AntiFam_results/
hmmsearch --cpu $THREADS --cut_ga -o "$DB_name"_hmmsearch_vs_"$input_name".tsv --domtblout domtblout_"$DB_name"_hmmsearch_vs_"$input_name".tsv --tblout tabular_"$DB_name"_hmmsearch_vs_"$input_name".tsv $DBpath $input_seqs
# Inspect results
