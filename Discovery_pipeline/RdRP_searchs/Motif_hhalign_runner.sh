#!/bin/bash
#hostname

##################################################################################################
usage() {
	echo ' 
      Arguments:
    #	Desc (suggestion) [default]
   -t	Threads (all) [2]
   -o	output directory (...) [pwd]
   -y input type (accepted "Multi" or "single", for alignments or single seq per sequence) [single]
   -i	Input (for MSAs, path to directory containing ONLY the MSAs (keep suffix as .afa or .faa). For multi fasta file, path to .faa file )
   -M	Motif type (#/ID) 
   -d	motifDB path [./Motif_DB/]
   -p	prefix for output .tsv file  [prefff]
   -s	suffix for output .tsv file  [sufff]
   -r	Remove temps? ([False]) 
   -e	E-value [(0.01)] 

   Usage: bash /home/neri/Documents/GitHub/RVMT/RdRp_Search/shell/Motif_hhalign_runner.sh -t 10  -o /media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/set_c210212/output/  -y "single"   -i /media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/set_c210212/singletons.faa -M RdRp304 -d /media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/set_c210212/hhdb/db  -p preff -s suff -r False -e 0.1
'
	exit
}
##################################################################################################

##### Set defaults: #####
THREADS=2
output_dir=$(pwd)
input_type='single'
input='input.faa'
Motif_type='Mot'
motifDB_path='Motif_DB'
out_pr='prefff'
out_su='sufff'
rm_tmp="False"
eval=0.001

##### Parse args: #####
while getopts t:o:e:y:i:p:s:r:M:d: flag; do
	case "${flag}" in
	t) THREADS=${OPTARG} ;;
	o) output_dir=${OPTARG} ;;
	y) input_type=${OPTARG} ;;
	i) input=${OPTARG} ;;
	M) Motif_type=${OPTARG} ;;
	d) motifDB_path=${OPTARG} ;;
	p) out_pr=${OPTARG} ;;
	s) out_su=${OPTARG} ;;
	r) rm_tmp=${OPTARG} ;;
	e) eval=${OPTARG} ;;
	*) usage ;;
	esac
done

# mandatory arguments
# if [ ! "$input" ] || [ ! "$VM" ]; then
if [ ! "$input" ]; then
	usage
	echo "arguments -i must be provided"
fi

mkdir $output_dir
cd $output_dir

if [ "$input_type" == "Multi" ]; then
	echo "Using alignments from $input"
	ini_name="$(echo "$input" | rev | cut -d"/" -f2 | rev)"
	mkdir HHMsearch_"$ini_name"_vs_mot"$Motif_type"
	input_fastas=$(echo $(ls "$input"))
	echo ${input_fastas} >tmparg_file.txt
	sed -i 's| |\n|g' tmparg_file.txt
	parallel -a tmparg_file.txt -j"$THREADS" hhsearch -i $input/{} -d $motifDB_path -o /dev/null -cpu 1 -e $eval -hide_cons -Z 100000 -show_ssconf -blasttab HHMsearch_"$ini_name"_vs_mot"$Motif_type"/$(basename {.} .afa)_hhsearch_rawout.tsv
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
	parallel -a tmparg_file.txt -j"$THREADS" hhsearch -i "$ini_name"_splitted/{} -d $motifDB_path -o /dev/null -cpu 1 -e $eval -hide_cons -Z 100000 -show_ssconf -blasttab HHMsearch_"$ini_name"_vs_mot"$Motif_type"/$(basename {.} .faa)_hhsearch_rawout.tsv
	cat HHMsearch_"$ini_name"_vs_mot"$Motif_type"/* >"$out_pr"_HHMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su"_rawout.tsv
fi

sed -i 's|.Cons.msa||g' "$out_pr"_HHMsearch_vs_"$ini_name"_vs_mot"$Motif_type"_"$out_su"_rawout.tsv # Comment this line out if you didn't use Neri/Wolf style suffix indicating hcon was used when creating the profiles.

if [ "$rm_tmp" != "False" ]; then
	rm HHMsearch_"$ini_name"_vs_mot"$Motif_type" "$ini_name"_splitted -r
fi
echo "Done"

echo "Printing params to $(pwd)/params.txt"

Nseqs=$(grep ">" $input -c)
echo "Motif_hhalign_runner.sh called on $Nseqs sequences at $(date) by: " >$(pwd)/params.txt

echo "$@" >>$(pwd)/params.txt
