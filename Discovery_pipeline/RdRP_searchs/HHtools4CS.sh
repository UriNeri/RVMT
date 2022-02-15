#!/bin/bash
#hostname

# hhsearch pdb70 kv7 genes.
module load hh-suite3
module load blast/blast-2.6.0
module load perl
export PATH=$PATH:/urigo/DBs/resources/
mmseqs2=/urigo/DBs/resources/mmseqs/bin/mmseqs

# THREADS=$1 #Number of threads.
# Memory=$2  # In GB
# set_ID=$3
# db_name=$4

Algo="hhblits"

round=2.03
V2_dir=/scratch200/urineri/V2/ # Path to main WD.
input="$V2_dir"/input/V2_6TranX_TNset.faa
OUTDIR="$V2_dir"/IRS/ #Path to output directory.
new_profiles=/scratch200/urineri/V2/profiles/set_"$set_ID"/
profileDir=$new_profiles/hhm_profiles/

pdb70=/urigo/DBs/soeding/hh-dbs/pdb70/pdb70 #_hhm_db
pfam_32=/urigo/DBs/soeding/hh-dbs/pfam_32/pfam
uniclust30=/urigo/DBs/soeding/hh-dbs/uniclust30/uniclust30_2018_08

if [ "$db_name" == "pdb70" ]; then
	DB_PATH="$pdb70"
fi
if [ "$db_name" == "pfam_32" ]; then
	DB_PATH="$pfam_32"
fi
if [ "$db_name" == "uniclust30" ]; then
	DB_PATH="$uniclust30"
fi

cd $OUTDIR
Currn_out="$OUTDIR"/"$Algo"_set_"$set_ID"_vs_"$db_name"
mkdir "$Currn_out"
cd $Currn_out
mkdir hhblits_out
cd hhblits_out
echo "Started $Algo for input: set_$set_ID"
for profilehhm in $profileDir/*.hhm; do
	profile_name=$(basename $profilehhm ".msa.faa.hhm")
	echo "Started $Algo for input: $profile_name"
	hhsearch -i $profilehhm -d $DB_PATH -blasttab ./hhblits_"$db_name"_vs_set_"$set_ID"_"$profile_name".tsv -o ./hhblits_"$db_name"_vs_set_"$set_ID"_"$profile_name".stdout -cpu $THREADS -maxmem $Memory
done
cat ./*.tsv >../hhblits_"$db_name"_vs_set_"$set_ID".tsv
sed -i 's|.msa.faa||g' ../hhblits_"$db_name"_vs_set_"$set_ID".tsv
rm ./*.tsv
for srchstdout in ./*.stdout; do
	grep "cluster" $srchstdout >greped.tsv
	sed -i '/^Query/d' ./greped.tsv
	sed -i '/^Command/d' ./greped.tsv
	cat ./greped.tsv >>../hhblits_"$db_name"_vs_set_"$set_ID"_stdouts.txt
done

echo "Done $Algo for set_$set_ID"
