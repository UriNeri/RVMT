#!/bin/bash
#hostname

##################################################################################################
if [[ $# -eq 0 ]]; then
	echo '   
   Arguments:
   #	Desc (suggestion)	
   1	Threads
   2	output directory
   3	Input fasta file (expected *.faa)
   4	Motif type (#) 
   5	motifDB path 
   6	prefix for output .tsv file 
   7	suffix for output .tsv file 
   8	Remove temp dir / MMseqs2DBs (False)
'
	exit
fi
##################################################################################################
THREADS=$1
output_dir=$2
input_fasta=$3
Motif_type=$4
motifDB_path=$5
out_pr=$6
out_su=$7
rm_tmp=$8

qcov=1
eval=0.5

ini_name=$(basename $input_fasta ".faa")
mkdir $output_dir
cd $output_dir
cp $input_fasta ./
input_fasta="$output_dir"/"$ini_name".faa
mkdir "$ini_name"_blastp_DB tmp
cd $ini_name"_blastp_DB"
makeblastdb -in $input_fasta -dbtype 'prot' -title $ini_name -out $ini_name
blastp_DB="$output_dir"/"$ini_name"_blastp_DB/"$ini_name"
cd $output_dir

cd tmp
for cluster in "$motifDB_path"/*Cons.*.faa; do # for cluster in "$motifDB_path"/*Cons.faa;
	cluster_name=$(basename $cluster ".msa.faa")
	echo "input cluster: $cluster_name "
	psiblast -word_size 2 -evalue $eval -max_target_seqs 100000 -threshold 9 -dbsize 20000000 -ignore_msa_master -in_msa $cluster -qcov_hsp_perc $qcov -num_threads $THREADS -out "$cluster_name"_pisblast_"$ini_name".tsv -db "$blastp_DB" -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore"
	# echo "Done $cluster_name vs $ini_name"
done
ls *_pisblast_*.tsv | xargs -I% sed 's/$/\t%/' % >../"$out_pr"_psiblast_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv # To concatenate the results (And add the profile to the output based on the filename)
cd ../
tbf=_pisblast_"$ini_name".tsv
sed -i "s|.Cons$tbf||g" "$out_pr"_psiblast_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv
# sed -i "s|.Cons.faa$tbf||g" "$out_pr"_psiblast_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv

if [ "$rm_tmp" == "True" ]; then
	rm tmp "$ini_name"_blastp_DB -rf
fi
