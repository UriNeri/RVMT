cdv2
cd absent*
iter=7
mkdir cated_absenties_iter"$iter"
cd cated_absenties_iter"$iter"
cp ../cated_absenties_iter"$iter".fasta ./
# the following splits the fasta into parts with equal number of sequences each. later, each of these portions-fastas is passed to it's own job. Total number of resulting files should be == total number of running jobs.
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%10770==0){file=sprintf("iter7%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' <cated_absenties_iter"$iter".fasta
rm cated_absenties_iter"$iter".fasta
j=1
for i in ./*.fasta; do
	mv $i ./cated_absenties_iter"$iter"_"$j".fasta
	j=$(($j + 1))
	echo $j
done
j=1
rpj=9 # RAM in GB per job
tpj=6 #Threads per job
for i in ./*.fasta; do
	qsub -v numys="$j",THREADS="$tpj" -l mem="$rpj"gb,nodes=1:ppn=$tpj -q gophna /urigo/urineri/scripts/RNA_Viruses/v2/shell_scripts/BlastN_nt_parl.sh
	j=$(($j + 1))
	echo $j
done
# cat BlastN_*.tsv > ../BlastN_cated_absenties_iter6_vs_nt_Pid75_qcov10_eval0.01.tsv

# qsub -v numys="$j" -lmem="$rpj"gb,pmem="$rpj"gb,vmem="$rpj"gb,pvmem="$rpj"gb,nodes=1:ppn=$tpj -q gophna /urigo/urineri/scripts/RNA_Viruses/v2/shell_scripts/BlastN_nt_parl.sh

# # #qdel'er (In case you fuck up).
# j=334199 # <number of first job>
# for i in ./*.fasta
# do
# qdel "$j" # Say bye bye!
# j=$(( $j + 1 ))
# echo $j
# done
