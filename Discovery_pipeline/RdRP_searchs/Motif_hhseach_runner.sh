#!/bin/bash
#hostname

##################################################################################################
if [[ $# -eq 0 ]]; then
	echo '   
   Arguments:
   #	Desc (suggestion)	
   1	Threads
   2	output directory
   3  input type (accepted "Multi" or "single", for alignments or single seq per sequence)
   4	Input (for MSAs, path to directory containing ONLY the MSAs (keep suffix as .afa or .faa (.hhm also works)). For multi fasta file, path to .faa file )
   5	Motif type (#) 
   6	motifDB path (db)
   7	prefix for output .tsv file 
   8	suffix for output .tsv file 
   9	Remove temp dir / MMseqs2DBs (False)
   10 extension (.hhm)
'
	exit
fi
##################################################################################################
THREADS=$1
output_dir=$2
input_type=$3
input=$4
Motif_type=$5
motifDB_path=$6
out_pr=$7
out_su=$8
rm_tmp=$9
extensionis=${10}
eval=0.5

mkdir $output_dir
cd $output_dir

if [ "$input_type" == "Multi" ]; then
	echo "Using alignments from $input"
	ini_name="$(echo "$input" | rev | cut -d"/" -f2 | rev)"
	mkdir HHMsearch_"$ini_name"_vs_mot"$Motif_type"
	input_fastas=$(echo $(ls "$input"))
	echo ${input_fastas} >tmparg_file.txt
	sed -i 's| |\n|g' tmparg_file.txt
	if [ "$extensionis" == ".hhm" ]; then
		parallel -a tmparg_file.txt -j"$THREADS" hhsearch -e $eval -i $input/{} -d $motifDB_path -o /dev/null -cpu 1 -hide_cons -Z 100000 -blasttab HHMsearch_"$ini_name"_vs_mot"$Motif_type"/$(basename {.} .afa)_hhsearch_rawout.tsv
	fi
	if [ "$extensionis" != ".hhm" ]; then
		parallel -a tmparg_file.txt -j"$THREADS" hhsearch -M 50 -e $eval -i $input/{} -d $motifDB_path -o /dev/null -cpu 1 -hide_cons -Z 100000 -blasttab HHMsearch_"$ini_name"_vs_mot"$Motif_type"/$(basename {.} .afa)_hhsearch_rawout.tsv
	fi
	#TO-DO: add file name to each line of each file
	cat HHMsearch_"$ini_name"_vs_mot"$Motif_type"/* >"$out_pr"_HHMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su"_rawout.tsv
fi

if [ "$input_type" == "single" ]; then
	echo "Splitting to seqs"
	ini_name=$(basename $input ".faa")
	cp $input ./
	input="$output_dir"/"$ini_name".faa
	mkdir "$ini_name"_splitted HHMsearch_"$ini_name"_vs_mot"$Motif_type"
	cd "$ini_name"_splitted
	splitfasta.pl $input -ext .faa
	cd ../
	input_fastas=$(echo $(ls "$ini_name"_splitted))
	echo ${input_fastas} >tmparg_file.txt
	sed -i 's| |\n|g' tmparg_file.txt
	parallel -a tmparg_file.txt -j"$THREADS" hhsearch -M 50 -i {} -d $motifDB_path -o /dev/null -cpu 1 -hide_cons -Z 100000 -blasttab HHMsearch_"$ini_name"_vs_mot"$Motif_type"/$(basename {.} .faa)_hhsearch_rawout.tsv
	#TO-DO: add file name to each line of each file
	cat HHMsearch_"$ini_name"_vs_mot"$Motif_type"/* >"$out_pr"_HHMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su"_rawout.tsv
fi

tbf=".Cons.msa"
sed -i "s|$tbf||g" "$out_pr"_HHMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su"_rawout.tsv

if [ "$rm_tmp" == "True" ]; then
	rm HHMsearch_"$ini_name"_vs_mot"$Motif_type" "$ini_name"_splitted -r
fi
echo "Done"
