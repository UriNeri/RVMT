#!/bin/bash
#hostname

##################################################################################################
if [[ $# -eq 0 ]]; then
	echo '   
   Arguments:
   #	Desc (suggestion)	
   1	Threads
   2	output directory
   3	Input (either fasta (expected *.faa) or MMseqs2DB path)
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

mkdir $output_dir
cd $output_dir

ini_name=$(basename $input)
extension="${ini_name##*.}"
if [[ "$extension" == "faa" ]]; then
	echo "$input is .faa, making foramtting as MMseqs2DB"
	ini_name=$(basename $input ".faa")
	cp $input ./
	input="$output_dir"/"$ini_name".faa
	mkdir "$ini_name"_MMseqs2DB "$ini_name"_vs_motif."$Motif_type" tmp
	mmseqs createdb $input "$ini_name"_MMseqs2DB/MMseqs2DB
	querydb="$ini_name"_MMseqs2DB/MMseqs2DB
else
	echo "$input isn't .faa, assuming MMseqs2DB "
	querydb="$input"
	ini_name=$(basename "$input")
	mkdir "$ini_name"_vs_motif."$Motif_type" tmp
fi

targetdb=$motifDB_path
resultsdb="$ini_name"_vs_motif."$Motif_type"/resultsdb
mmseqs search $querydb $targetdb $resultsdb ./tmp/ -k 6 -s 7.5 -e $eval --threads $THREADS --split-memory-limit 100G --max-seqs 6000
mmseqs convertalis $querydb $targetdb $resultsdb "$out_pr"_MMseqs2_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su".tsv --format-output "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov"
if [ "$rm_tmp" == "True" ]; then
	rm tmp "$ini_name"_vs_motif."$Motif_type" -r
	if [[ "$extension" == "faa" ]]; then
		rm tmp "$ini_name"_MMseqs2DB/ -rf
	fi
fi
echo "Done"
