#!/bin/bash
#hostname

# module load blast/blast-2.6.0
# module load R/R-3.3.2
# module load gcc/gcc-7.2.0 python/python-anaconda3.5
export PATH=$PATH:/urigo/urineri/resources/new_bbmap/bbmap/
export PATH=$PATH:/urigo/urineri/scripts/RNA_Viruses/Neri_scripts/resources/
####### TO DO: #######
# 1. Add print check variable types.
# 2. Print if input has exception to log.err (by >> as it may already exist).
# DONE: 3. Check CPU type - if instruction avx2 exists opt to it, if not - report and switch to sse2.
# 4.

####### enviroment #######
THREADS=48
MEMORY=350g
Neri_scripts=/urigo/urineri/scripts/RNA_Viruses/Neri_scripts/
resources_dir="$Neri_scripts"/resources
V2_dir=/scratch200/urineri/V2/
# V2_dir=/davidb/urineri/V2/

export PATH=$PATH:"$resources_dir"/mmseqs_AVX2/bin/
# export PATH=$PATH:"$resources_dir"/mmseqs_SSE4/bin/
diamond=/urigo/DBs/resources/Vadaimond
Diamond_nr=/urigo/DBs/Diamond_nr/nr.dmnd
taxonmap=/urigo/DBs/nodamp/prot.accession2taxid.gz
BLASTn_NT=/urigo/DBs/nt/nt
step=8
workdir="$V2_dir"/MTs/step_"$step"/
refdirs="$V2_dir"/MGs/

# outdir="$workdir"/results_dir/

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

############ Filtering bits ############
# study_list=("s2" "s4")
# for i in ${study_list[*]};
# do
# cd $refdirs/"$i"
# mkdir MMseqs2DB
# echo "Started formating MMseqs2DB for $i"
# mmseqs createdb ./cated_sk100_reps.fna ./MMseqs2DB/MMseqs2DB # --dbtype 0 # --compressed 1
# # mmseqs  createindex ./MMseqs2DB/MMseqs2DB ./tmp  --threads $THREADS  --search-type 3
# echo "Finished formating MMseqs2DB for $i"
# # cd $MGs
# done
############ VS DNAome #################

# s_list=("s1" "s2" "s3")
# s_list=("s2") #Start easy - set checkpoint A.
# s="s4"
# for s in ${s_list[*]}
# do
# echo "Started working on $s "
# targetdb=$refdirs"/"$s"/"MMseqs2DB/MMseqs2DB
# for i in {2..3};
# do
# cd $workdir
# echo "Started working on C$i"
# cd C"$i"

# mkdir C"$i"_vs_"$s" C"$i"DB tmp
# # mmseqs createdb ./*.fasta ./C"$i"DB/C"$i"DB    #--dbtype 0 #--compressed 1
# querydb=./C"$i"DB/C"$i"DB
# resultsdb=C"$i"_vs_"$s"/resultsdb
# mmseqs search $querydb $targetdb $resultsdb ./tmp/   --threads $THREADS --search-type 3 -s 1  --min-seq-id 0.70  --min-aln-len 100 -e 0.000001 --max-seqs 50  --max-accept 1 --split-memory-limit 374G # --kmer-per-seq-scale 0.4  --compressed 1 # --cov-mode  -c 0.5 --start-sens 1 --sens-steps 1 -s 4
# # mmseqs search $querydb $targetdb $resultsdb ./tmp/ --force-reuse --threads $THREADS --search-type 3 -s 3  --min-seq-id 0.70  --min-aln-len 100 -e 0.0001 --max-seqs 75  --max-accept 1 --split-memory-limit 374G # --kmer-per-seq-scale 0.4  --compressed 1 # --cov-mode  -c 0.5 --start-sens 1 --sens-steps 1 -s 4
# mmseqs convertalis $querydb $targetdb $resultsdb MMseqs2_C"$i"_vs_"$s".tsv --search-type 3 --format-output "query,target,evalue,pident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,mismatch"
# cd $workdir
# echo "Done with C$i"
# done
# echo "Done with $s"
# done
# Filter hits in postprocessing - resplit the filtered set to (fewer) 14.3GB independent chunks, and reformat MMseqs2 DBs for these.

# ############# VS IMG/VR (DNA only) ################
# s="VR_DNA"
# echo "Started filtering using $s "
# targetdb="$V2_dir"IMG_VR/MMseqs2DB/MMseqs2DB
# cd $workdir
# echo "Started working on $s"
# mkdir step_2_vs_"$s"
# querydb=./step_2DB/step_2DB
# resultsdb=step_2_vs_"$s"/resultsdb
# mmseqs search $querydb $targetdb $resultsdb ./tmp/ --threads $THREADS --search-type 3 -s 2  --min-seq-id 0.70  --min-aln-len 100 -e 0.00001 --max-seqs 50  --max-accept 1 # --kmer-per-seq-scale 0.4  --compressed 1 # --cov-mode  -c 0.5 --start-sens 1 --sens-steps 1 -s 4
# mmseqs convertalis $querydb $targetdb $resultsdb MMseqs2_step_2_vs_"$s".tsv --search-type 3 --format-output "query,target,evalue,pident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,mismatch"
# echo "Done filtering with $s"
# # Filter hits in postprocessing - resplit the filtered set to (fewer) 14.3GB independent chunks, and reformat DBs for these.

# ############ VS NCBI's NR  ############
# s="NR"
# Algo="DiamondX"
# pid=75
# evalue=0.00001
# echo "Started $Algo using $s on evalue = $evalue ; pid = $pid % "
# cd $workdir
# input_seq=DNAVR_fitered.fasta
# mkdir tmp

# echo "Started filtering input: $(basename $input_seq .fasta)"
# $diamond blastx -b5.0 -p $THREADS -k 2 -d $Diamond_nr -q $input_seq --evalue $evalue -o DiamondX_vs_"$s".tsv --tmpdir ./tmp/ --outfmt 6 qseqid sseqid staxids stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles --taxonmap $taxonmap -v
# echo "Done filtering with $s"
# # Filter hits in postprocessing - resplit the filtered set to (fewer) 14.3GB independent chunks, and reformat DBs for these.

# # ############ VS FPs DB  ############
# s="FPs"
# Algo="MMseqs2"
# cd $workdir
# mkdir tmp step_"$step"_vs_"$s"
# targetdb="$V2_dir"/FPs/MMseqs2_DB/MMseqsDB
# querydb=./MMseqs2DB/MMseqs2DB
# resultsdb=step_"$step"_vs_"$s"/resultsdb
# mmseqs search $querydb $targetdb $resultsdb ./tmp/   --threads $THREADS --search-type 3 -s 1  --min-seq-id 0.70  --min-aln-len 120 -e 0.000001 --max-seqs 50  --max-accept 1 --split-memory-limit 374G # --kmer-per-seq-scale 0.4  --compressed 1 # --cov-mode  -c 0.5 --start-sens 1 --sens-steps 1 -s 4
# # # mmseqs search $querydb $targetdb $resultsdb ./tmp/ --force-reuse --threads $THREADS --search-type 3 -s 3  --min-seq-id 0.70  --min-aln-len 100 -e 0.0001 --max-seqs 75  --max-accept 1 --split-memory-limit 374G # --kmer-per-seq-scale 0.4  --compressed 1 # --cov-mode  -c 0.5 --start-sens 1 --sens-steps 1 -s 4
# mmseqs convertalis $querydb $targetdb $resultsdb MMseqs2_vs_"$s".tsv --search-type 3 --format-output "query,target,evalue,pident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,mismatch"
# # cd $workdir
# echo "Done filtering with $s"
# Filter hits in postprocessing - resplit the filtered set to (fewer) 14.3GB independent chunks, and reformat DBs for these.

# # # ############ VS Head S2 (s5) DB  ############
# s="1"
# step=6
# workdir="$V2_dir"/MTs/step_"$step"/
# Algo="MMseqs2"
# cd $workdir
# mkdir tmp step_"$step"_vs_s"$s"
# targetdb="$V2_dir"/MGs/s1/MMseqs2DB/MMseqs2DB
# querydb=./MMseqs2DB/MMseqs2DB
# resultsdb=step_"$step"_vs_s"$s"/resultsdb
# mmseqs search $querydb $targetdb $resultsdb ./tmp/   --threads $THREADS --force-reuse --search-type 2 -s 2  --min-seq-id 0.9  --min-aln-len 33 -e 0.000001 --max-seqs 50  --max-accept 1 --split-memory-limit 374G # --kmer-per-seq-scale 0.4  --compressed 1 # --cov-mode  -c 0.5 --start-sens 1 --sens-steps 1 -s 4
# # # mmseqs search $querydb $targetdb $resultsdb ./tmp/ --force-reuse --threads $THREADS --search-type 3 -s 3  --min-seq-id 0.70  --min-aln-len 100 -e 0.0001 --max-seqs 75  --max-accept 1 --split-memory-limit 374G # --kmer-per-seq-scale 0.4  --compressed 1 # --cov-mode  -c 0.5 --start-sens 1 --sens-steps 1 -s 4
# mmseqs convertalis $querydb $targetdb $resultsdb MMseqs2_vs_"$s".tsv --search-type 2 --format-output "query,target,evalue,pident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,mismatch"
# # cd $workdir
# echo "Done filtering with $s"
# Filter hits in postprocessing - resplit the filtered set to (fewer) 14.3GB independent chunks, and reformat DBs for these.

# # ############ VS NCBI's NT  ############
# s="NT"
# cd $workdir
# Algo="BLASTn"
# pid=65
# evalue=0.000001
# echo "Started $Algo using $s on evalue = $evalue ; pid = $pid % "
# export BLASTDB=$BLASTDB:/urigo/DBs/nt
# input_seq=DiamondX_MG_filtered.fasta
# echo "Started filtering input: $(basename $input_seq .fasta)"
# blastn -task megablast -perc_identity $pid -num_threads $THREADS -max_target_seqs 1 -query $input_seq -db "$BLASTn_NT" -evalue $evalue -out BLASTn_vs_"$s".tsv -outfmt "6 qseqid sseqid sscinames pident qlen length mismatch gapope evalue bitscore stitle sskingdoms staxids"
# cd $workdir
# echo "Done filtering with $s"
# # # Filter hits in postprocessing - resplit the filtered set to (fewer) 14.3GB independent chunks, and reformat DBs for these.

# # # ############ mmseqs2 VS Pfam32_A  ############
# s="Pfam32_A"
# Algo="MMseqs2"
# cd $workdir
# mkdir tmp step_"$step"_vs_"$s"
# targetdb="/urigo/DBs/soeding/mmseq2_dbs/pfam/pfam_profile"
# resultsdb=./step_"$step"_vs_"$s"/resultsdb
# querydb="$workdir"/MMseqs2DB/MMseqs2DB

# mmseqs search $querydb $targetdb $resultsdb tmp -k 5 -s 5 -e 0.001  --threads $THREADS # --force-reuse
# mmseqs convertalis $querydb $targetdb $resultsdb result_tab  --format-output "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov"
# # # # rm ./querydb ./tmp ./resultsdb -r
# # # module load R/R-3.3.2
# # # Rscript /urigo/DBs/soeding/mmseq2_dbs/mmseq_pfam32_parser.r $(pwd)

# # # # ############ mmseqs2 VS FPs  ############
# # s="FPDB"
# # Algo="MMseqs2"
# # cd $workdir
# # mkdir tmp step_"$step"_vs_"$s"
# # targetdb="new_FPDB/FPDB"
# # resultsdb=./resultdb/resultdb
# # querydb=./repeatfiltered1DB/MMseqs2DB
# # mmseqs search $querydb $targetdb $resultsdb tmp  -s 3  --min-seq-id 0.6  --min-aln-len 100 -e 0.0000001  --threads 10 --search-type 3  --max-seqs 50  --max-accept 1 # --force-reuse
# # mmseqs convertalis $querydb $targetdb $resultsdb vsrepeatFPdb.tsv --search-type 3 --format-output "query,target,evalue,pident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,mismatch"
# # # # rm ./querydb ./tmp ./resultsdb -r
# # # module load R/R-3.3.2
# # # Rscript /urigo/DBs/soeding/mmseq2_dbs/mmseq_pfam32_parser.r $(pwd)

# mmseqs linsearch $querydb $targetdb $resultsdb ./tmp/ --threads $THREADS --search-type 3   --min-seq-id 0.75  --min-aln-len 119 -e 0.0000001 --max-accept 1 --split-memory-limit 100G # --kmer-per-seq-scale 0.4  --compressed 1 # --cov-mode  -c 0.5 --start-sens 1 --sens-steps 1 -s 4
# mmseqs easy-linsearch sk100_rep_seq.fasta  $refdirs/s1/cated_sk100_reps.fna $resultsdb ./tmp/ --threads $THREADS --search-type 3   --min-seq-id 0.75  --min-aln-len 119 -e 0.0000001 --max-accept 1 --split-memory-limit 100G

# ############ Filtering bits ############
# # study_list=("s1" "s2" "s3" "s4" "s5" "s6")
# scov=70
# qcov=50
# pid=-1
# evalue=0.00001

# input_seq="$workdir"/DiamondX_MG_filtered.fasta
# study_list=("s1" "s2" "s3" "s4")

# for i in ${study_list[*]};
# do
# cd $workdir
# targetdb=$refdirs/"$i"/dmnd/"$i".dmnd
# echo "Started diamond blastx vs $i"
# # diamond makedb --in cated_sk100_reps.faa  --threads $THREADS  --db dmnd/"$i".dmnd
# diamond blastx -b4.0  -p $THREADS -k 1 --query-cover $qcov  --evalue $evalue -q $input_seq -d $targetdb -o step_"$step"_vs_"$i".tsv -f 6 qseqid sseqid qstart qend sstart send evalue bitscore length pident mismatch slen #--subject-cover $scov
# # callgenes.sh in=cated_sk100_reps.fna outa=cated_sk100_reps.faa  threads=$THREADS  overwrite=true #-Xmx"$MEMORY"
# # mmseqs createdb ./cated_sk100_reps.fna ./MMseqs2DB/MMseqs2DB # --dbtype 0 # --compressed 1
# # mmseqs  createindex ./MMseqs2DB/MMseqs2DB ./tmp  --threads $THREADS  --search-type 3
# echo "Finished diamond blastx vs $i"
# # cd $MGs
# done

cd "$V2_dir"/MTs/
gunzip -k V3_contigs.fasta.gz
