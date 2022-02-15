#!/bin/bash
#hostname

##################################################################################################
if [[ $# -eq 0 ]]; then
	echo '
   Create 6 frame end to end translation of a given nucelic fasta file. Output written to where the input was.
   Arguments:
   1    input fasta file
   2    Output to stout (Default = False). Only  "True" makes the output go to stdout.
'
	exit
fi
##################################################################################################
# export PATH=$PATH:/urigo/urineri/resources/perl/
export PATH=$PATH:/home/neri/Documents/GitHub/seals-2/bin/ali/
export PATH=$PATH:/home/neri/Documents/GitHub/seals-2/bin/convert/
export PATH=$PATH:/home/neri/Documents/GitHub/seals-2/bin/stat/
export PATH=$PATH:/home/neri/Documents/GitHub/seals-2/bin/tab/
export PATH=$PATH:/home/neri/Documents/GitHub/seals-2/bin/tree/
export PATH=$PATH:/home/neri/Documents/GitHub/seals-2/bin/misc/
export PATH=$PATH:/home/neri/Documents/GitHub/seals-2/YIW/
# module load perl
# cd /scratch200/urineri/V3/input/
# input="/scratch200/urineri/V3/V301/input/4set2/scafs4set2search.fasta"
STOUD_is=$2
input=$1

if [ $STOUD_is == "True" ]; then
	cat "$input" | orf -c=1 -f=0 -m=0 -n=2 -w=1 | fa_strict -l=1 -r=X -u
elif [ $STOUD_is != "True" ]; then
	input_name="$(basename "$input")"
	input_name="$(echo "$input_name" | rev | cut -f 2- -d '.' | rev)"
	dir_name=$(dirname "$input")
	cat "$input" | orf -c=1 -f=0 -m=0 -n=2 -w=1 | fa_strict -l=1 -r=X -u >"$dir_name"/"$input_name"_sixframe.faa

fi
# cat $input |  fa_select -p="^\d+" -w=0 |  orf -c=1 -f=0 -m=0 -n=2 -w=1 | perl -pe 's/^(>\S+)\.(-*\d+)/\1\_\2/' |  fa_strict -l=1 -r=X -u > "$input_name"_sixframe.faa
# cat $input |  fa_select -p="^\d+" -w=0 |  orf -c=1 -f=0 -m=0 -n=2 -w=1  |  fa_strict -l=1 -r=X -u > "$dir_name"/"$input_name"_sixframe.faa
# cat "$input" |  orf -c=1 -f=0 -m=0 -n=2 -w=1  |  fa_strict -l=1 -r=X -u > "$dir_name"/"$input_name"_sixframe.faa
# cat RdRp_scafs.fasta |  fa_select -p="^\d+" -w=0 |  orf -c=1 -f=0 -m=0 -n=2 -w=1 | perl -pe 's/^(>\S+)\.(-*\d+)/\1\_\2/' |  fa_strict -l=1 -r=X -u > V3_sixframe.fa
