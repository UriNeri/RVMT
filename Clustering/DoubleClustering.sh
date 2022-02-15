#!/bin/bash
#hostname

THREADS=10
work_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre/"
# input_RdRps="$work_dir""VR0520.faa"
input_RdRps="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre/DoubleClustering/ALL_RdRp_raw.faa"
ini_name=$(basename $input_RdRps ".faa")
cd $work_dir || exit
mkdir DoubleClustering && cd DoubleClustering  || exit

#### First Clustering at 0.9id ####
mkdir "$ini_name"_MMseqs2DB_1 DB_clu_rep_50 DB_clu_rep_98 MM85clusters MM50clusters 
querydb="$ini_name"_MMseqs2DB_1/MMseqs2DB
mmseqs createdb "$input_RdRps" "$querydb"
mmseqs cluster "$querydb" MM85clusters/clusterDB tmp -c 0.333 -e 0.1 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.9  --threads $THREADS 
mmseqs result2repseq "$querydb" MM85clusters/clusterDB DB_clu_rep_98/DB_clu_rep_98
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_98/DB_clu_rep_98 clu_reps98.faa  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM85clusters/clusterDB results98DB_clu.tsv
mmseqs easy-linclust $input_RdRps sk98.clu tmp  --min-seq-id 0.9 -c 0.9 --cov-mode 1 --kmer-per-seq-scale 0.4 --threads $THREADS 
#### Second Clustering at 0.5id ####
input_RdRps="clu_reps98.faa"
ini_name=$(basename $input_RdRps ".faa")
mkdir "$ini_name"_MMseqs2DB_2/
querydb="$ini_name"_MMseqs2DB_2/MMseqs2DB
mmseqs createdb $input_RdRps "$querydb"
querydb="$ini_name"_MMseqs2DB_2/MMseqs2DB
mmseqs cluster "$querydb" MM50clusters/clusterDB tmp -c 0.333 -e 0.1 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.5  --threads $THREADS 
mmseqs result2repseq "$querydb" MM50clusters/clusterDB DB_clu_rep_50/DB_clu_rep_50
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_50/DB_clu_rep_50 clu_reps50.faa  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM50clusters/clusterDB resultsDB_clu.tsv


#### Unrelated Clustering at 0.98id ####
THREADS=10
work_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/Clsts2/"
input_Contigs="$work_dir""IDFT.fasta"
ini_name=$(basename $input_Contigs ".fasta")
cd $work_dir 

mkdir "$ini_name"_MMseqs2DB_1  DB_clu_rep_98 MM85clusters  
querydb="$ini_name"_MMseqs2DB_1/MMseqs2DB
mmseqs createdb "$input_Contigs" "$querydb"
mmseqs cluster "$querydb" MM85clusters/clusterDB tmp -c 0.333 -e 0.001 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.9  --threads $THREADS 
mmseqs result2repseq "$querydb" MM85clusters/clusterDB DB_clu_rep_98/DB_clu_rep_98
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_98/DB_clu_rep_98 clu_reps98.fasta  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM85clusters/clusterDB resultsDB_clu.tsv

mkdir "$ini_name"  DB_clu_rep_80 MM80clusters  
# querydb="$ini_name"_MMseqs2DB_1/MMseqs2DB
# mmseqs createdb "$input_Contigs" "$querydb"
mmseqs cluster "$querydb" MM80clusters/clusterDB tmp -c 0.333 -e 0.001 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.8  --threads $THREADS 
mmseqs result2repseq "$querydb" MM80clusters/clusterDB DB_clu_rep_80/DB_clu_rep_80
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_80/DB_clu_rep_80 clu_reps80.fasta  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM80clusters/clusterDB resultsDB_clu80.tsv


#####
precent_id=0.9 # <0.99 AND >0.95

mmseqs easy-linclust $ini_name".fasta" sk98.clu tmp --min-seq-id 0.9 -c 0.9 --cluster-mode 2 --cov-mode 0 --kmer-per-seq-scale 0.4 --threads $THREADS 
mmseqs easy-linclust $ini_name".fasta" sk100.clu tmp --min-seq-id 0.99 -c 1.0 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS 
mmseqs easy-linclust $ini_name".fasta" sk99.clu tmp  --min-seq-id 0.9 -c 0.9 --cov-mode 1 --kmer-per-seq-scale 0.4 --threads $THREADS 




#### Nth Clustering at 0.9id and 0.5 ####
THREADS=10
work_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre_3/Double_Clustering3/"
input_RdRps="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre_3/Double_Clustering3/Hybrid_RdRps.faa"
ini_name=$(basename $input_RdRps ".faa")
cd $work_dir 
mkdir DoubleClustering && cd DoubleClustering  
mkdir "$ini_name"_MMseqs2DB_1 DB_clu_rep_50 DB_clu_rep_98 MM85clusters MM50clusters 
querydb="$ini_name"_MMseqs2DB_1/MMseqs2DB
mmseqs createdb "$input_RdRps" "$querydb"
mmseqs cluster "$querydb" MM85clusters/clusterDB tmp -c 0.333 -e 0.1 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.9  --threads $THREADS 
mmseqs result2repseq "$querydb" MM85clusters/clusterDB DB_clu_rep_98/DB_clu_rep_98
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_98/DB_clu_rep_98 clu_reps98.faa  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM85clusters/clusterDB results98DB_clu.tsv
mmseqs easy-linclust $input_RdRps sk98.clu tmp  --min-seq-id 0.9 -c 0.9 --cov-mode 1 --kmer-per-seq-scale 0.4 --threads $THREADS 

#### Second Clustering at 0.5id ####
input_RdRps="clu_reps98.faa"
ini_name=$(basename $input_RdRps ".faa")
mkdir "$ini_name"_MMseqs2DB_2/
querydb="$ini_name"_MMseqs2DB_2/MMseqs2DB
mmseqs createdb $input_RdRps "$querydb"
querydb="$ini_name"_MMseqs2DB_2/MMseqs2DB
mmseqs cluster "$querydb" MM50clusters/clusterDB tmp -c 0.333 -e 0.1 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.5  --threads $THREADS 
mmseqs result2repseq "$querydb" MM50clusters/clusterDB DB_clu_rep_50/DB_clu_rep_50
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_50/DB_clu_rep_50 clu_reps50.faa  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM50clusters/clusterDB resultsDB_clu.tsv


#### Nth+1 Clustering at 0.85id and 0.5 ####
THREADS=10
work_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre_3/Double_Clustering3/"
input_RdRps="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre_3/Double_Clustering3/Raw_And_Moded_RdRp.faa"
ini_name=$(basename $input_RdRps ".faa")
cd $work_dir
mkdir DoubleClustering && cd DoubleClustering 
mkdir "$ini_name"_MMseqs2DB_1 DB_clu_rep_50 DB_clu_rep_85 MM85clusters MM50clusters 
querydb="$ini_name"_MMseqs2DB_1/MMseqs2DB
mmseqs createdb "$input_RdRps" "$querydb"
mmseqs cluster "$querydb" MM85clusters/clusterDB tmp -c 0.333 -e 0.1 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.85  --threads $THREADS 
mmseqs result2repseq "$querydb" MM85clusters/clusterDB DB_clu_rep_85/DB_clu_rep_85
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_85/DB_clu_rep_85 clu_reps85.faa  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM85clusters/clusterDB results85DB_clu.tsv
mmseqs easy-linclust $input_RdRps sk85.clu tmp  --min-seq-id 0.85 -c 0.85 --cov-mode 1 --kmer-per-seq-scale 0.4 --threads $THREADS 
#### Second Clustering at 0.5id ####
input_RdRps="clu_reps85.faa"
ini_name=$(basename $input_RdRps ".faa")
mkdir "$ini_name"_MMseqs2DB_2/
querydb="$ini_name"_MMseqs2DB_2/MMseqs2DB
mmseqs createdb $input_RdRps "$querydb"
querydb="$ini_name"_MMseqs2DB_2/MMseqs2DB
mmseqs cluster "$querydb" MM50clusters/clusterDB tmp -c 0.333 -e 0.1 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.5  --threads $THREADS 
mmseqs result2repseq "$querydb" MM50clusters/clusterDB DB_clu_rep_50/DB_clu_rep_50
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_50/DB_clu_rep_50 clu_reps50.faa  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM50clusters/clusterDB resultsDB_clu.tsv




#### 10.10.2021 - Unrelated Clustering at 0.98id ####
THREADS=12
work_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre_3/Double_Clustering3/Nuc98/"
input_Contigs="/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/ALL_nuc_3007.fasta"
ini_name=$(basename $input_Contigs ".fasta")
cd $work_dir 

mkdir "$ini_name"_MMseqs2DB_1  DB_clu_rep_98 MM95clusters  
querydb="$ini_name"_MMseqs2DB_1/MMseqs2DB
mmseqs createdb "$input_Contigs" "$querydb"
mmseqs cluster "$querydb" MM95clusters/clusterDB tmp -c 0.75 -e 0.0000000001 --cov-mode 5 --cluster-mode 2 --min-seq-id 0.98  --gap-open 7 --gap-extend 3 --threads $THREADS 
mmseqs result2repseq "$querydb" MM95clusters/clusterDB DB_clu_rep_98/DB_clu_rep_98
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_98/DB_clu_rep_98 clu_reps98.fasta  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM95clusters/clusterDB resultsDB_clu.tsv

mmseqs easy-linclust $ini_name".fasta" sk98.clu tmp --gap-open 7 --gap-extend 3 --min-seq-id 0.98 -e 0.0000000001 -c 1.0 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS 

mmseqs easy-linclust $input_Contigs sk99.clu tmp --min-seq-id 0.99 -e 0.0000000001 -c 0.99 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS 

mmseqs easy-linclust $input_Contigs sk99.clu tmp --min-seq-id 0.99 -e 0.000000000001 -c 0.99 --cov-mode 1   --threads $THREADS 

mmseqs easy-linclust $input_Contigs sk99.clu tmp --min-seq-id 0.95 -e 0.00000000000000001 -c 1 --cluster-mode 2 --cov-mode 1   --threads $THREADS --seq-id-mode 2



#### Nth Clustering at 0.9id and 0.5 ####
THREADS=10
work_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre_3/Double_Clustering4/"
mkdir $work_dir
input_RdRps="/media/HDD1/uri/RNA_Vir_MTs/V3/Vre_3/Double_Clustering3/Hybrid_RdRps2.faa"
ini_name=$(basename $input_RdRps ".faa")
cd $work_dir 
mkdir DoubleClustering && cd DoubleClustering  
mkdir "$ini_name"_MMseqs2DB_1 DB_clu_rep_50 DB_clu_rep_90 MM90clusters MM50clusters 
querydb="$ini_name"_MMseqs2DB_1/MMseqs2DB
mmseqs createdb "$input_RdRps" "$querydb"
mmseqs cluster "$querydb" MM90clusters/clusterDB tmp -c 0.333 -e 0.0000000001 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.9  --threads $THREADS 
mmseqs result2repseq "$querydb" MM90clusters/clusterDB DB_clu_rep_90/DB_clu_rep_90
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_90/DB_clu_rep_90 clu_reps90.faa  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM90clusters/clusterDB results90DB_clu.tsv
mmseqs easy-linclust $input_RdRps sk90.clu tmp  --min-seq-id 0.9 -c 0.9 --cov-mode 1 --kmer-per-seq-scale 0.4 --threads $THREADS 
#### Second Clustering at 0.5id ####
input_RdRps="clu_reps90.faa"
ini_name=$(basename $input_RdRps ".faa")
mkdir "$ini_name"_MMseqs2DB_2/
querydb="$ini_name"_MMseqs2DB_2/MMseqs2DB
mmseqs createdb $input_RdRps "$querydb"
querydb="$ini_name"_MMseqs2DB_2/MMseqs2DB
mmseqs cluster "$querydb" MM50clusters/clusterDB tmp -c 0.333 -e 0.1 --cov-mode 1 --cluster-mode 2 --min-seq-id 0.5  --threads $THREADS 
mmseqs result2repseq "$querydb" MM50clusters/clusterDB DB_clu_rep_50/DB_clu_rep_50
mmseqs result2flat "$querydb" "$querydb" DB_clu_rep_50/DB_clu_rep_50 clu_reps50.faa  --use-fasta-header
mmseqs createtsv "$querydb" "$querydb" MM50clust
blastn -task dc-megablast  -evalue 0.001 -perc_identity 66 -word_size 12  -max_target_seqs 10000000 -num_threads $THREADS -db V302filt_DB -out /urigo/urineri/DCblastn_IDFT_vs_V302filt.tsv   -outfmt "6 qseqid sseqid pident sstart send qstart qend slen qlen length evalue bitscore"
ers/clusterDB resultsDB_clu.tsv

mmseqs easy-linclust $input_RdRps sk90.clu tmp --min-seq-id 0.90 -e 0.00000000001 -c 0.90 --cluster-mode 2 --cov-mode 1   --threads $THREADS --seq-id-mode 2
mmseqs easy-linclust $input_RdRps sk50.clu tmp --min-seq-id 0.50 -e 0.00000000001 -c 0.75 --cluster-mode 2 --cov-mode 1   --threads $THREADS --seq-id-mode 2
mmseqs easy-linclust $input_RdRps sk75.clu tmp --min-seq-id 0.75 -e 0.00000000001 -c 0.85 --cluster-mode 2 --cov-mode 1   --threads $THREADS --seq-id-mode 2
mmseqs easy-linclust $input_RdRps sk85.clu tmp --min-seq-id 0.85 -e 0.00000000001 -c 0.85 --cluster-mode 2 --cov-mode 1   --threads $THREADS --seq-id-mode 2
mmseqs easy-linclust $input_RdRps sk80.clu tmp --min-seq-id 0.80 -e 0.00000000001 -c 0.75 --cluster-mode 2 --cov-mode 1   --threads $THREADS --seq-id-mode 2



input=/media/HDD1/uri/RNA_Vir_MTs/V3/Simon/IDFT101021.fna

mmseqs easy-linclust $input sk98.clu tmp --min-seq-id 0.98 -e 0.00000000001 -c 0.98 --cluster-mode 2 --cov-mode 1   --threads $THREADS --seq-id-mode 2



input=/media/HDD1/uri/RNA_Vir_MTs/V3/Simon/All.101021.fna
THREADS=11
mmseqs easy-linclust $input sk90.clu tmp  --min-seq-id 0.90 -e 0.00000000001 -c 0.90 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS 
mmseqs easy-linclust $input NC95.clu tmp  --min-seq-id 0.95 -e 0.00001 -c 0.95 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS 
mmseqs easy-linclust $input NC98.clu tmp --min-seq-id 0.98 -e 0.00001 -c 0.98 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS 



# For tANI 08.10.2021
cd /media/HDD1/uri/RNA_Vir_MTs/V3/Antonio/tANI/AVA/

input=/media/HDD1/uri/RNA_Vir_MTs/V3/Simon/All.101021.fna
THREADS=11
mmseqs easy-search $input $input All101021.AVA.blout  tmp  --threads $THREADS -s 7.5 --search-type 3



THREADS=48
input=/scratch200/urineri/AvA/IDFT.fna
mmseqs easy-search $input $input IDFT.AVA.blout   tmp  --add-self-matches true --threads $THREADS  -s 5.7 --search-type 3 -e 0.00001 --disk-space-limit 1T --max-seqs 5000000 --min-aln-len 100 --format-output query,target,fident,nident,pident,qstart,qend,tstart,tend,tlen,qlen,alnlen,mismatch,gapopen,evalue,bits,tcov,qcov

input=/scratch200/urineri/AvA/All.101021.fna
mmseqs easy-search $input $input All101021.AVA.blout  tmp  --add-self-matches true --threads $THREADS  -s 5.7 --search-type 3 -e 0.00001 --disk-space-limit 1T --max-seqs 5000000 --min-aln-len 100 --format-output query,target,fident,nident,pident,qstart,qend,tstart,tend,tlen,qlen,alnlen,mismatch,gapopen,evalue,bits,tcov,qcov



 awk -F "\t" '{ if(($12 >= 18) && ($5 >= 40 && $18  >= 0.2 && $17  >= 0.2)) { print } }' All101021.AVA.blout  > Awked_All101021.AVA.blout 



c("qid", "tid", "fid", "nid", "pid", "q1", "q2", "t1", "t2", "tL", "qL", "alnlen", "mismatch", "gapopen", "evalue", "bits", "tcov", "qcov")
  1       2       3     4       5     6      7    8     9     10    11     12          13         14          15     16       17      18

 pv All101021.AVA.blout | awk -F "\t" '{ if(($16 >= 20) && ($5 >= 60)) { print } }' > Awked_All101021.AVA.blout 
1,2,5,10,11,12,15,16,17,18
pv All101021.AVA.blout | awk -F "\t" '{ if(($16 >= 20) && ($5 >= 60)) { print } }' | cut -f1-2,4-12,15- > Awked_All101021.AVA.blout 

c("qid", "tid", "fident", "pident", "q1", "q2", "t1", "t2", "tL", "qL", "alnlen", "mismatch", "gapopen", "evalue", "bits", "tcov", "qcov")

THREADS=48
input=/scratch200/urineri/AvA/All.101021.fna
mmseqs easy-search $input $input All101021.AVA.blout  tmp --min-seq-id 0.5 -a true --add-self-matches true --threads $THREADS  -s 3 --search-type 3 -e 0.0000001 --disk-space-limit 1T --max-seqs 5000000 --min-aln-len 100 --format-output query,target,nident,pident,qstart,qend,tstart,tend,tlen,qlen,alnlen,mismatch,gapopen,evalue,bits,tcov,qcov


      qid           tid          fid    nid     pid    q1   q2    t1    t2      tL       qL     alnlen    mismatch   gapopen     evalue     bits    tcov    qcov
ND_0525700      ND_0525700      1.000   838     100.   1    838    1    838     838     838     838           0        0       0.000E+00    1503    1.000   1.000


c("qid", "tid", "nid", "pid", "q1", "q2", "t1", "t2", "tL", "qL", "alnlen",  "evalue", "bits", "tcov", "qcov")
  1        2       3     4     5     6      7    8     9     10    11         12          13         14          15     16       17      18
 pv All101021.AVA.blout | awk -F "\t" '{ if(($16 >= 20) && ($5 >= 60)) { print } }' > Awked_All101021.AVA.blout 
 pv Awked_All101021.AVA.blout | awk -F "\t" '{ if((($12 + 0) < 1E-6) && ($4 >= 75)) { print } }' > DoubleAwked_All101021.AVA.blout


 pv Awked_All101021.AVA.blout | awk -F "\t" '{ if((($12 + 0) < 1E-6) && ($4 >= 75)) { print } }' > DoubleAwked_All101021.AVA.blout



c("qid", "tid", "fid", "nid", "pid", "q1", "q2", "t1", "t2", "tL", "qL", "alnlen", "mismatch", "gapopen", "evalue", "bits", "tcov", "qcov")
  1       2       3     4       5     6      7    8     9     10    11     12          13         14          15     16       17      18

pv All101021.AVA.blout | awk -F "\t" '{ if(($16 >= 20) && (($15 + 0) < 1E-6) && ($5 >= 89)) { print } }' | cut -f1-2,4-12,15- > Awked_All101021.AVA.blout 
)



THREADS=48
mmseqs easy-search ../MissingNDs.fna ../All.101021.fna MissingNDs.vs.All.101021.blout  tmp --min-seq-id 0.5 -a true --add-self-matches true --threads $THREADS  -s 3 --search-type 3 -e 0.0000001 --disk-space-limit 1T --max-seqs 5000000 --min-aln-len 100 --format-output query,target,nident,pident,qstart,qend,tstart,tend,tlen,qlen,alnlen,evalue,bits,tcov,qcov
mmseqs easy-search ../All.101021.fna ../MissingNDs.fna All.101021.vs.MissingNDs.blout  tmp --min-seq-id 0.5 -a true --add-self-matches true --threads $THREADS  -s 3 --search-type 3 -e 0.0000001 --disk-space-limit 1T --max-seqs 5000000 --min-aln-len 100 --format-output query,target,nident,pident,qstart,qend,tstart,tend,tlen,qlen,alnlen,evalue,bits,tcov,qcov
mmseqs easy-search ../MissingNDs.fna ../MissingNDs.fna MissingNDs.vs.MissingNDs.blout  tmp --min-seq-id 0.5 -a true --add-self-matches true --threads $THREADS  -s 3 --search-type 3 -e 0.0000001 --disk-space-limit 1T --max-seqs 5000000 --min-aln-len 100 --format-output query,target,nident,pident,qstart,qend,tstart,tend,tlen,qlen,alnlen,evalue,bits,tcov,qcov


pv All101021.AVA.blout | awk -F "\t" '{ if(($16 >= 20) && (($15 + 0) < 1E-6) && ($5 >= 89)) { print } }' | cut -f1-2,4-12,15- > Awked_All101021.AVA.blout 



c("qid", "tid",  "nid", "pid", "q1", "q2", "t1", "t2", "tL", "qL", "alnlen",  "evalue", "bits", "tcov", "qcov")
  1       2       3     4       5     6      7    8     9     10     11          12          13         14          15     16       17      18

pv MissingNDs.blout | awk -F "\t" '{ if(($13 >= 20) && (($12 + 0) < 1E-6) && ($4 >= 89)) { print } }' MissingNDs.blout > Awked_MissingNDs.blout 

mmseqs easy-search ../MissingNDs.fna ../MissingNDs.fna MissingNDs.vs.MissingNDs.blout  tmp --min-seq-id 0.5 -a true --add-self-matches true --threads $THREADS  -s 3 --search-type 3 -e 0.0000001 --disk-space-limit 1T --max-seqs 5000000 --min-aln-len 100 --format-output query,target,nident,pident,qstart,qend,tstart,tend,tlen,qlen,alnlen,evalue,bits,tcov,qcov
