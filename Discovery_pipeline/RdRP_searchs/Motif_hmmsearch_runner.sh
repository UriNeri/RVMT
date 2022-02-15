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
input=$3
Motif_type=$4
motifDB_path=$5
out_pr=$6
out_su=$7
rm_tmp=$8
eval=0.5

ini_name=$(basename $input ".faa")
echo $ini_name
mkdir $output_dir
cd $output_dir
cp $input ./
input="$output_dir"/"$ini_name".faa
mkdir "$ini_name"_hmmsearch_motif."$Motif_type"/
cd "$ini_name"_hmmsearch_motif."$Motif_type"/
motifDB_path_name=$(basename $motifDB_path ".hmm") #removes the path of the RdRP profile.
echo "input motifDB_path: $motifDB_path_name "
hmmsearch --noali --cpu $THREADS -E $eval --incE $eval -o "$out_pr"_HMMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su"_rawout.tsv --domtblout "$out_pr"_HMMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv --tblout tabular_"$out_pr"_HMMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv $motifDB_path $input
sed -i '/^#/d' "$out_pr"_HMMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv
echo "Done $motifDB_path_name vs $ini_name"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"}' "$out_pr"_HMMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv >../"$out_pr"_HMMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv #Legacy - don't chage to OFS because of trailing sep issue with parser
cd ../

if [ "$rm_tmp" == "True" ]; then
	rm tmp "$ini_name"_hmmsearch_motif."$Motif_type"/ -r
fi
echo "Done"
