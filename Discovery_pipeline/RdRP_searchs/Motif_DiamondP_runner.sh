#!/bin/bash
#hostname

##################################################################################################
if [[ $# -eq 0 ]]; then
	echo '   
   Arguments:
   #	Desc (suggestion)	
   1	Threads
   2	output directory
   3	Input (expected *.faa)
   4	Motif type (#) 
   5	motifDB path (singleton.faa)
   6	prefix for output .tsv file 
   7	suffix for output .tsv file 
   8	Remove temp dir / dmnd DB (False)
'
	exit
fi
##################################################################################################
THREADS=$1
output_dir=$2
input=$3
Motif_type=$4
motifDB_path=$5
out_pr=$6
out_su=$7
rm_tmp=$8
scov=10
qcov=50
eval=0.05
min_id=60
mkdir $output_dir
cd $output_dir
ini_name=$(basename $input)
echo "$input is .faa, making foramtting as dmnd DB"
ini_name=$(basename $input ".faa")
cp $input ./
input="$output_dir"/"$ini_name".faa
mkdir "$ini_name"_Diamond_DB
diamond makedb --in $input -d "$ini_name"_Diamond_DB/"$ini_name".dmnd -p $THREADS
Diamond_DB="$ini_name"_Diamond_DB/"$ini_name".dmnd
diamond blastp -q $motifDB_path --db $Diamond_DB -b3.0 -k 0 -p $THREADS --query-cover $qcov --id $min_id --more-sensitive --evalue $eval -o "$out_pr"_DiamondP_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv --outfmt 6 sseqid qseqid evalue gapopen pident length qstart qend sstart send qlen slen mismatch bitscore -v
awk -F'\t' -vOFS='\t' -v Motif_type="$Motif_type" '{ $2 = "mot."Motif_type"."$2}1' "$out_pr"_DiamondP_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv >tmpout
mv tmpout "$out_pr"_DiamondP_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv
if [ "$rm_tmp" == "True" ]; then
	rm tmp "$ini_name"_Diamond_DB -r
fi
echo "Done"
