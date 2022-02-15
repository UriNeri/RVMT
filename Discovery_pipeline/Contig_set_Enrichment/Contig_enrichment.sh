#!/bin/bash
#hostname

export PATH=$PATH:"/urigo/DBs/resources/tools/mmseqs-current/bin/"
THREADS=48
mem=340g

Min_Ali=120
Min_ID=0.66
Sens=1
evalue=0.000000001

tmpDir=/scratch200/urineri/tmp/

cd /scratch200/urineri/ReSearch/

mmseqs createdb --dbtype 2 SimonDFT.fasta IDFTdb/IDFT
queryDB=/scratch200/urineri/ReSearch/IDFTdb/IDFT

mkdir StudiesDB
mkdir output
outdb=/scratch200/urineri/ReSearch/output/outdb

mmseqs createdb --dbtype 2 Cated_MTs.fasta StudiesDB/Studies
targetDB=/scratch200/urineri/ReSearch/StudiesDB/Studies

mmseqs search --search-type 3 -e $evalue --min-aln-len $Min_Ali --min-seq-id $Min_ID -s $Sens -c 0.85 --cov-mode 1 --max-seqs 1000000 --remove-tmp-files TRUE -a TRUE --threads $THREADS --disk-space-limit 1T --split-memory-limit 330G $queryDB $targetDB $outdb $tmpDir

mmseqs convertalis --search-type 3 $queryDB $targetDB $outdb /scratch200/urineri/ReSearch/SimonDFT_vsMTs.tsv --format-output "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,mismatch,qcov,tcov"
