#!/bin/bash
#hostname

export PATH=$PATH:/urigo/DBs/resources/tools/bbmap/
export PATH=$PATH:/urigo/DBs/resources/bin/
alias diamond='/urigo/DBs/resources/bin/diamond'

export PATH=$PATH:/urigo/DBs/resources/tools/hmmer-current/
export PATH=$PATH:/urigo/DBs/resources/tools/hmmer-current/bin/
export PATH=$PATH:/urigo/DBs/resources/tools/hhsuite3.3/bin/
export PATH=$PATH:/urigo/DBs/resources/tools/hhsuite3.3/
export PATH=$PATH:/urigo/DBs/resources/tools/mmseqs-current/
export PATH=$PATH:/urigo/DBs/resources/tools/mmseqs-current/bin/
# export  PATH=$PATH:"/urigo/DBs/resources/tools/ncbi-current/bin"

for file in alimask hmmalign hmmbuild hmmconvert hmmemit hmmfetch hmmlogo hmmpgmd hmmpgmd_shard hmmpress hmmscan hmmsearch hmmsim hmmstat jackhmmer phmmer nhmmer nhmmscan makehmmerdb; do
	tpb=$(echo /urigo/DBs/resources/tools/hmmer-current/bin/$file)
	alias $file=$(echo $tpb)
done

####### Global enviroment #######
THREADS=48       # uncomment if not passed via the qsub (-v) command
MEMORY=340g      # uncomment if not passed via the qsub (-v) command
algo="hmmsearch" # uncomment if not passed via the qsub (-v) command (hhsuite sepretaly)
eval=0.1
scov=1
qcov=1

input=/scratch200/urineri/Annotations/All3007.faa
input_dir=$(dirname $input)
ini_name=$(basename $input ".faa")
OUTDIR=/scratch200/urineri/Annotations/Search_Out/hmmsearch_vs_CDD3.19/
mkdir $OUTDIR
profileDir=/scratch200/urineri/lys/C11

if [ "$algo" == "hmmsearch" ]; then
	echo "Running $algo"
	###### hmmsearch (HMMER) ######
	# module load hmmr/hmmr-3.1b2
	hmmdb=$profileDir/*.hmm # HMM format (HMMER3 !!!! not hh-suite).

	mkdir $OUTDIR/hmmsearch
	cd $OUTDIR/hmmsearch

	hmmdb_name=$(basename $hmmdb ".hmm")
	echo "input hmmdb: $hmmdb_name "
	# hmmsearch  --noali   --cpu $THREADS -E $eval --incE $eval -o "$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --domtblout domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --tblout tabular_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv $hmmdb $input
	# hmmsearch  --noali   --cpu $THREADS -E $eval --incE $eval -o "$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --domtblout domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --tblout tabular_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv $hmmdb $input
	eval=0.01

	# hmmsearch  --max --cpu $THREADS -E $eval --domE $eval --incdomE $eval --incE $eval -o "$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --domtblout domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --tblout tabular_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv $hmmdb $input
	hmmsearch --cpu $THREADS -E $eval --domE $eval --incdomE $eval --incE $eval -o "$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --domtblout domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --tblout tabular_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv $hmmdb $input

	sed -i '/^#/d' domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv
	echo "target_name	target_accession	qL	query_name	query_accession	pL	E-value_fullseq	score_fullseq	bias_fullseq	#	of	c-Evalue	i-Evalue	bias_thisdomain	score_thisdomain	p1	p2	q1	q2	env_from	env_to	acc	description_of_target" >awked_domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv >>awked_domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv #Legacy - don't chage to OFS because of trailing sep issue with parser
	echo "Done $hmmdb_name vs  $ini_name"
fi

if [ "$algo" == "psiblast" ]; then
	echo "Running $algo"

	##### PSI-BLAST (blast-2.6.0) #####
	module load blast/blast-2.6.0

	cd "$input_dir"/
	mkdir $ini_name"_blastp_DB"
	cd $ini_name"_blastp_DB"
	makeblastdb -in $input -dbtype 'prot' -title $ini_name -out $ini_name
	blastp_DB="$input_dir"/"$ini_name"_blastp_DB/"$ini_name"

	mkdir $OUTDIR/psiblast_vs_set_"$ini_name"
	cd $OUTDIR/psiblast_vs_set_"$ini_name"
	cluster=/scratch200/urineri/lys/C4/msaFiles/RapidLysis_awked.afa
	for cluster in $profileDir/MSA/*.fas; do
		cluster_name=$(basename $cluster ".fas")
		echo "input cluster: $cluster_name "
		# psiblast -word_size 2 -evalue $eval -max_target_seqs 10000000 -threshold 9 -dbsize 20000000 -ignore_msa_master -in_msa $cluster -qcov_hsp_perc $qcov -num_threads $THREADS -out "$cluster_name"_pisblast_vs_"$ini_name".tsv -db "$blastp_DB" -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore"
		psiblast -word_size 2 -dbsize 22000000 -max_target_seqs 10000000 -threshold 9 -in_msa $cluster -num_threads $THREADS -out "$cluster_name"_pisblast_vs_"$ini_name".tsv -db "$blastp_DB" -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore"

		echo "Done $cluster_name vs $ini_name"
	done
	ls *_pisblast_vs*.tsv | xargs -I% sed 's/$/\t%/' % >../cated_psiblast_vs_V"$ini_name".tsv # To concatenate the results (And add the profile to the output based on the filename)
	cd ../
	tbf=.Cons_pisblast_vs_"$ini_name".tsv
	sed -i "s/$tbf//g" cated_psiblast_vs_V"$ini_name".tsv
	tbf=_pisblast_vs_"$ini_name".tsv
	sed -i "s/$tbf//g" cated_psiblast_vs_V"$ini_name".tsv

	# rm psiblast_vs_set_"$set_ID" -r
	echo "Done PSI-BLAST vs  $ini_name"
fi

if [ "$algo" == "mmseqs2" ]; then
	echo "Running $algo"
	##### MMseqs2 (profile search) ######
	# module load mmseq2
	rm /scratch200/urineri/tmp/ -r
	cd "$input_dir"
	mkdir "$ini_name"_MMseqs2DB
	mmseqs createdb $input "$ini_name"_MMseqs2DB/MMseqs2DB
	querydb="$input_dir"/"$ini_name"_MMseqs2DB/MMseqs2DB

	cd $OUTDIR
	mkdir "$ini_name"_profiles_vs_MMseqs2DB tmp

	seqdb="$profileDir"/MMseqs2_profiles/Cls.Lys_MM
	resultsdb="$ini_name"_profiles_vs_MMseqs2DB/resultsdb
	mmseqs search $querydb $seqdb $resultsdb /scratch200/urineri/tmp/ -s 7.5 -e $eval --threads $THREADS
	mmseqs convertalis $querydb $seqdb $resultsdb "$ini_name"_profiles_vs_MMseqs2DB.tsv --format-output "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov"
	rm /scratch200/urineri/tmp/ "$ini_name"_profiles_vs_MMseqs2DB -r
fi

if [ "$algo" == "diamond" ]; then
	echo "Running $algo"
	###### Diamond BLASTp (Diamond-Ver xxx) ######
	rm /scratch200/urineri/tmp1/ -r
	mkdir /scratch200/urineri/tmp1/

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

	mkdir $OUTDIR/DiamondP_vs_set_"$ini_name"
	cd $OUTDIR/DiamondP_vs_set_"$ini_name"
	singletons="$profileDir"/MegaMerged.fa
	echo "Started DiamondP $ini_name"
	diamond blastp --ultra-sensitive -q $singletons --db $Diamond_DB -c1 -b5.0 -k 1000000000 --culling-overlap 1 -p $THREADS -t /scratch200/urineri/tmp1/ --evalue $eval -o ./DiamondP_singltons_vs_"$ini_name".tsv --outfmt 6 qseqid sseqid pident length qstart qend sstart send qlen slen mismatch gapopen evalue bitscore -v
	diamond blastp --ultra-sensitive --db $singletons -q $input -c1 -b5.0 -k 1000000000 --culling-overlap 1 -p $THREADS -t /scratch200/urineri/tmp1/ --evalue $eval -o ./DiamondP_"$ini_name"_vs_singltons.tsv --outfmt 6 qseqid sseqid pident length qstart qend sstart send qlen slen mismatch gapopen evalue bitscore -v

	# --outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore"
	echo "Done DiamondP $ini_name"
	rm /scratch200/urineri/tmp1/ -r
fi

if [ "$algo" == "jackhmmer" ]; then
	echo "Running $algo"
	###### jackhmmer (HMMER) ######
	seqdb="$profileDir"/MegaMerged.fa
	eval=0.00001
	iter=2
	mkdir $OUTDIR/jackhmmer
	cd $OUTDIR/jackhmmer

	seqdb_name=$(basename $seqdb ".fa")
	echo "input hmmdb: $seqdb_name "
	# hmmsearch  --noali   --cpu $THREADS -E $eval --incE $eval -o "$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --domtblout domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --tblout tabular_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv $hmmdb $input
	# hmmsearch  --noali   --cpu $THREADS -E $eval --incE $eval -o "$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --domtblout domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --tblout tabular_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv $hmmdb $input
	# hmmsearch  --max --cpu $THREADS -E $eval --incE $eval -o "$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --domtblout domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv --tblout tabular_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv $hmmdb $input
	jackhmmer -N $iter --incE $eval --cpu $THREADS -A "$seqdb_name"_jackhmmer_vs_"$ini_name".ali -o "$seqdb_name"_jackhmmer_vs_"$ini_name".tsv --domtblout domtblout_"$seqdb_name"_jackhmmer_"$ini_name".tsv --tblout tabular_"$seqdb_name"_jackhmmer_"$ini_name".tsv $seqdb $input

	# sed -i '/^#/d' domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv
	# awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"}' domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv > awked_domtblout_"$hmmdb_name"_hmmsearch_vs_"$ini_name".tsv #Legacy - don't chage to OFS because of trailing sep issue with parser
	echo "Done $seqdb_name vs  $ini_name"
fi
