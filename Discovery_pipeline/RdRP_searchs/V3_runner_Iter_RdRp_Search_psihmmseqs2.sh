#!/bin/bash
#hostname

export PATH=$PATH:/urigo/urineri/resources/new_bbmap/bbmap/
export PATH=$PATH:/urigo/urineri/scripts/RNA_Viruses/Neri_scripts/resources/
####### TO DO: #######
# 2. Print if input has exception to log.err (by >> as it may already exist).

####### Global enviroment #######
THREADS=48     # uncomment if not passed via the qsub (-v) command
MEMORY=350g    # uncomment if not passed via the qsub (-v) command
algo=hmmsearch # uncomment if not passed via the qsub (-v) command (hhsuite sepretaly)
Neri_scripts=/urigo/urineri/scripts/RNA_Viruses/Neri_scripts/
resources_dir="$Neri_scripts"/../resources
V3_dir=/scratch200/urineri/V3/

export PATH=$PATH:"$resources_dir"/mmseqs_AVX2/bin/

# avx2_count=$(cat /proc/cpuinfo | grep avx2 -c)
# if [ "$avx2_count" == "0" ]
# then
#  print "No avx2 instruction - opting to sse2/4"
#  export PATH=$PATH:"$resources_dir"/mmseqs_SSE4/bin/
# fi
# if [ "$avx2_count" != "0" ] # prefered !
# then
# print "positive for avx2 instruction - opting to avx2"
# export PATH=$PATH:"$resources_dir"/mmseqs_AVX2/bin/
# fi

round=simon
set_ID=2
# V3_dir=/davidb/urineri/V3/
V3_dir=/scratch200/urineri/V3/V"$round"
input_dir="$V3_dir"/input/4set0
input="$input_dir"/for_set2_sixframe.faa
# input="$input_dir"/fungi_transcripts_nm_sixframe.faa
# input="$input_dir"/vr0520.faa
OUTDIR="$V3_dir"/output/set_"$set_ID" #Path to output directory.

new_profiles="$V3_dir"/profiles/set_"$set_ID"/
profileDir=$new_profiles #Legacy, change.
ini_name=$(basename $input ".faa")

if [ "$algo" == "hmmsearch" ]; then
	echo "Running $algo"
	###### hmmsearch (HMMER) ######
	module load hmmr/hmmr-3.1b2
	eval=0.5

	hmmdb=$profileDir/set_"$set_ID"_profiles.hmm # HMM format (HMMER3 engine !!!! not hh-suite).
	cd $OUTDIR
	mkdir hmmsearch
	hmmdb_name=$(basename $hmmdb | cut -d'.' -f 1) #removes the path of the RdRP profile.
	echo "input hmmdb: $hmmdb_name "
	# hmmsearch  --noali --max  --cpu $THREADS -E 0.5 --incE 0.5 -o $OUTDIR/"$hmmdb_name"_hmmsearch_vs_V"$round".tsv --domtblout $OUTDIR/domtblout_"$hmmdb_name"_hmmsearch_vs_V"$round".tsv --tblout $OUTDIR/tabular_"$hmmdb_name"_hmmsearch_vs_V"$round".tsv $hmmdb $input
	hmmsearch --noali --cpu $THREADS -E $eval --incE $eval -o $OUTDIR/hmmsearch/"$hmmdb_name"_hmmsearch_vs_V"$round".tsv --domtblout $OUTDIR/hmmsearch/domtblout_"$hmmdb_name"_hmmsearch_vs_V"$round".tsv --tblout $OUTDIR/hmmsearch/tabular_"$hmmdb_name"_hmmsearch_vs_V"$round".tsv $hmmdb $input
	cd hmmsearch
	sed -i '/^#/d' domtblout_set_"$set_ID"_profiles_hmmsearch*
	echo "Done $hmmdb_name vs V$round  && set_$set_ID"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"}' domtblout_set_"$set_ID"_profiles_hmmsearch_vs_V"$round".tsv >awked_domtblout_set_"$set_ID"_profiles_hmmsearch_vs_V"$round".tsv #Legacy - don't chage to OFS because of trailing sep issue with parser
fi

if [ "$algo" == "psiblast" ]; then
	echo "Running $algo"
	##### PSI-BLAST (blast-2.6.0) #####
	module load blast/blast-2.6.0
	qcov=5
	eval=0.5

	cd "$input_dir"/
	mkdir $ini_name"_blastp_DB"
	cd $ini_name"_blastp_DB"
	makeblastdb -in $input -dbtype 'prot' -title $ini_name -out $ini_name
	blastp_DB="$input_dir"/"$ini_name"_blastp_DB/"$ini_name"

	cd $OUTDIR
	mkdir psiblast_vs_set_"$set_ID"
	cd psiblast_vs_set_"$set_ID"
	for cluster in $profileDir/msaFiles/*.faa; do
		cluster_name=$(basename $cluster ".msa.faa")
		echo "input cluster: $cluster_name "
		psiblast -word_size 2 -evalue $eval -max_target_seqs 10000000 -threshold 9 -dbsize 20000000 -ignore_msa_master -in_msa $cluster -qcov_hsp_perc $qcov -num_threads $THREADS -out "$cluster_name"_pisblast_V"$round"_set_"$set_ID".tsv -db "$blastp_DB" -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore"
		echo "Done $cluster_name vs V$round  && set_$set_ID"
	done
	ls *_pisblast_V*.tsv | xargs -I% sed 's/$/\t%/' % >../cated_psiblast_vs_V"$round"_set_"$set_ID".tsv # To concatenate the results (And add the profile to the output based on the filename)
	cd ../
	tbf=_pisblast_V"$round"_set_"$set_ID".tsv
	sed -i "s/$tbf//g" cated_psiblast_vs_V"$round"_set_"$set_ID".tsv

	# rm psiblast_vs_set_"$set_ID" -r
	echo "Done PSI-BLAST vs V$round  && set_$set_ID"
fi

if [ "$algo" == "mmseqs2" ]; then
	echo "Running $algo"
	##### MMseqs2 (profile search) ######
	module load hh-suite3
	# module load mmseq2
	rm /scratch200/urineri/tmp/ -r
	eval=0.5
	cd "$input_dir"
	mkdir "$ini_name"_MMseqs2DB
	mmseqs createdb $input "$ini_name"_MMseqs2DB/MMseqs2DB
	querydb="$input_dir"/"$ini_name"_MMseqs2DB/MMseqs2DB

	cd $OUTDIR
	mkdir set_"$set_ID"_profiles_vs_V3 tmp

	targetdb="$profileDir"/MMseqs2_profiles/MMseqs2DB
	resultsdb=set_"$set_ID"_profiles_vs_V3/resultsdb
	mmseqs search $querydb $targetdb $resultsdb /scratch200/urineri/tmp/ -k 6 -s 7.5 -e $eval --threads $THREADS
	mmseqs convertalis $querydb $targetdb $resultsdb MMseqs2_vs_V"$round"_set_"$set_ID".tsv --format-output "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov"
	rm /scratch200/urineri/tmp/ "$OUTDIR"/set_"$set_ID"_profiles_vs_V3/ -r
fi

if [ "$algo" == "diamond" ]; then
	echo "Running $algo"
	###### Diamond BLASTp (Diamond-Ver xxx) ######
	rm /scratch200/urineri/tmp1/ -r
	mkdir /scratch200/urineri/tmp1/
	module load diamond
	scov=10
	qcov=10
	eval=0.5
	# if test -f "$Diamond_DB";
	# then
	#     echo "$Diamond_DB exists";
	# else
	# echo "$Diamond_DB does not exist, proceding to make it";
	# mkdir $Diamond_DB;
	# diamond makedb --in $input --db "$Diamond_DB"/"$ini_name".dmnd;
	# fi

	cd "$input_dir"
	mkdir "$ini_name"_Diamond_DB
	diamond makedb --in $input -d "$ini_name"_Diamond_DB/"$ini_name".dmnd -p $THREADS
	Diamond_DB="$input_dir"/"$ini_name"_Diamond_DB/"$ini_name".dmnd

	cd $OUTDIR
	mkdir DiamondP_vs_set_"$set_ID"
	cd DiamondP_vs_set_"$set_ID"
	singletons="$profileDir"/singletons.faa
	echo "Started DiamondP set_""$set_ID" " singletons vs V$round"
	# diamond blastp -q $singleton --db $Diamond_DB -b3.0 -p $THREADS  --subject-cover $scov --evalue $eval -o DiamondP_set_"$set_ID"_singltons_vs_V"$round".tsv --threads $THREADS --outfmt 6 -v
	diamond blastp -q $singletons --db $Diamond_DB -b4.0 -k 0 -p $THREADS -t /scratch200/urineri/tmp1/ --query-cover $qcov --more-sensitive --evalue $eval -o ./DiamondP_set_"$set_ID"_singltons_vs_V"$round".tsv --outfmt 6 qseqid sseqid pident length qstart qend sstart send qlen slen mismatch gapopen evalue bitscore -v
	# --outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore"
	echo "Done DiamondP set_""$set_ID" " singletons vs V$round"
	rm /scratch200/urineri/tmp1/ -r
fi
