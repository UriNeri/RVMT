#!/bin/bash
#hostname

THREADS=48

Neri_scripts=/urigo/urineri/scripts/RNA_Viruses/Neri_scripts/
resources_dir="$Neri_scripts"/resources
V2_dir=/scratch200/urineri/V2/
diamond=/urigo/DBs/resources/Vadaimond
export PATH=$PATH:"$resources_dir"/mmseqs_AVX2/bin/
Diamond_nr=/urigo/DBs/Diamond_nr/nr.dmnd
taxonmap=/urigo/DBs/nodamp/prot.accession2taxid.gz
BLASTn_NT=/urigo/DBs/nt/nt
step=2
workdir="$V2_dir"/MTs/step_"$step"/
refdirs="$V2_dir"/MGs/
# outdir="$workdir"/results_dir/

# avx2_count=$(cat /proc/cpuinfo | grep avx2 -c)
# if [ "$avx2_count" == "0" ]
# then
#  print "No avx2 instruction - opting to sse2"
#  export PATH=$PATH:"$resources_dir"/mmseqs_SSE4/bin/
# fi
# if [ "$avx2_count" != "0" ] # prefered !
# then
# print "positive for avx2 instruction - opting to avx2"
# export PATH=$PATH:"$resources_dir"/mmseqs_AVX2/bin/
# fi

############ mkdirs and DBs #################
for i in {0..3}; do
	cd $workdir
	echo "Started working on C$i"
	cd C"$i"
	mkdir ./C"$i"DB
	querydb=./C"$i"DB/C"$i"DB
	# resultsdb=C"$i"_vs_"$s"/resultsdb
	mmseqs createdb ./step_1_filtered_ml1200."$i".fasta $querydb #--dbtype 0 --compressed 1  #--split-memory-limit 32G
	cd $workdir
	echo "Done with C$i"
done
