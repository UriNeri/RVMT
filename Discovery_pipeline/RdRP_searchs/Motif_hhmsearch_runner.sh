HHM_runner="/media/uri/HDD1/uri/RNA_Vir_MTs/Neri_scripts/RdRp_search/shell/Motif_HHalign_runner.sh"

THREADS=10
output_dir=/media/uri/HDD1/uri/RNA_Vir_MTs/V3/V301/Motifs/Iterative/Iterations/Inflation_1.6/Iter_0/mot.1/
input_type="Multi"
input=/media/uri/HDD1/uri/RNA_Vir_MTs/V3/V301/ali.200525.05/ali.200525.05/
Motif_type=1
motifDB_path=/media/uri/HDD1/uri/RNA_Vir_MTs/V3/V301/Motifs/Iterative/Iterations/Inflation_1.6/Iter_0/mot.1/hhdb/db
out_pr="prfx"
out_su="sfx"
rm_tmp=True
eval=0.5

# cd /media/uri/HDD1/uri/RNA_Vir_MTs/V3/V301/ali.200525.05/
# for i in ./ali.200525.05/*.afa
# do
# file_with_suffix=$(basename "$i")
# file_name=$(basename $file_with_suffix ".afa")
# cat $i | fa2sr -w=0 |sr_filter -conplus -hcon=0 -ncon="$file_name"_con | sr2fa > ./tmphhm/"$file_name".Cons.msa.fas
# # hmmbuild --informat afa -n $file_name -o logs/"$file_name".hmm.log HMMfiles/"$file_name".hmm  ./msaFiles/"$file_name".Cons.msa.faa
# hhmake -v 2 -name $file_name -i ./tmphhm/"$file_name".Cons.msa.fas -o ./tmphhm/"$file_name".hhm -v 0  -id 100 -diff inf -M first
# done
