#!/usr/bin/env bash
module load hh-suite3
module load mmseq2
# THREADS=40
round=2.03

# Wolf="/urigo/urineri/RNA_Viruses/Wolf/new/"
# RdRp_DB="$Wolf"/MMseqs2_profiles/profileDb

# V2_dir="/urigo/urineri/RNA_Viruses/V2/results_V2_3/"
# V2_seqs="$V2_dir"V2_Scaffolds.fasta
# cd $V2_dir
# mkdir RdPs_vs_V2 MMseqs2_V2_DB
# # mmseqs createdb $V2_seqs "$V2_dir"/MMseqs2_DB/V2_results
# mmseqs createdb $V2_seqs "$V2_dir"/MMseqs2_V2_DB/V2_DB

# targetdb="$RdRp_DB"
# querydb="$V2_dir"/MMseqs2_V2_DB/V2_DB
# resultsdb="RdPs_vs_V2/results"
# mmseqs search $querydb $targetdb $resultsdb /scratch200/urineri/tmp/ -k 6 -s 7.5 -e 0.01  --threads $THREADS
# mmseqs convertalis $querydb $targetdb $resultsdb MMseqs2_search_results.tsv  --format-output "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov"
# rm  /scratch200/urineri/tmp/ RdPs_vs_V2/ "$V2_dir"/MMseqs2_V2_DB/ -r

# ###### MMseqs2 vs nt (for taxonomy report)  ######
THREADS=30
echo "Started MMseqs2 vs nt "
V2_dir="/urigo/urineri/RNA_Viruses/V2/results_V2_3/"
V2_New_hits_seqs="$V2_dir"V2_scaffolds.fasta
workdir=/scratch200/urineri/
cd $workdir
mkdir New_hits_nr MMseqs2_V2_New_hits_seqs_DB

querydb="$workdir"/MMseqs2_V2_New_hits_seqs_DB/New_hits_DB
mmseqs createdb $V2_New_hits_seqs "$querydb"

resultsdb="New_hits_nr/resultsdb
"targetdb="/scratch200/urineri/nr/nr.faaDB"
mmseqs search $querydb $targetdb $resultsdb /scratch200/urineri/tmp/ --threads $THREADS -s 1
mmseqs convertalis $querydb $targetdb $resultsdb MMseqs2_nr_vs_V"$round".tsv --format-output "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov"
mmseqs createtaxdb nr/nr.faaDB nr.faa.taxidmapping /urigo/DBs/NCBI_nr/taxonomy/ /scratch200/urineri/tmp/

# mmseqs createtsv $querydb $resultsdb taxonomyResult.tsv
# mmseqs taxonomyreport $targetdb $resultsdb taxonomyResult_report.ks # .ks = Kraken-style
# mmseqs taxonomyreport $targetdb $resultsdb Krona_taxonomy_report.html --report-mode 1
# rm ./querydb ./tmp ./taxonomyResult -r
echo "Finshed taxonomy w/ report"
