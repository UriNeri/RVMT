cd /media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/ORFs

seq="/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/ALL_nuc_3007.fasta"
prodigal -n -a ./trans_file_anon_mode_n.faa -d ./nuc_file_anon_mode_n.fna -i $seq -p meta -o ./verb_anon_mode_n.gff -f "gff"

grep "seqhdr=" ./verb_anon_mode.gff >seq_hdrs.txt
grep "transl_table=" ./verb_anon_mode.gff >ttbls.txt
awk -F ';' '{print $3}' seq_hdrs.txt | sed 's|seqhdr=||g' | sed 's|"||g' >seq_hdrs2.txt
awk -F ';' '{print $5}' ttbls.txt | sed 's|transl_table=||g' >ttbls2.txt

paste seq_hdrs2.txt ttbls2.txt >Gene_tables_anon_node.txt
rm ttbls2.txt seq_hdrs2.txt seq_hdrs.txt ttbls.txt

seq="/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/ORFs/ALL_nuc_3007_t4.fasta"
prodigal -i $seq -g 4 -o ./verb_t4.gff -f "gff" -a ./trans_file_t4.faa -d ./nuc_file_at4.fna

cd /media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/workdir2

Code_List=(6 14 16 22)
for i in ${Code_List[*]}; do # Loops over a directory of individual SRA runs as seperate subdirs.
	seq="/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/workdir2/IDFT_"$i".fasta"
	prodigal -i $seq -g $i -o ./verb_t_"$i".gff -f "gff" -a ./trans_file_t_"$i".faa -d ./nuc_file_t_"$i".fna
done

Code_List=(16 22)
for i in ${Code_List[*]}; do # Loops over a directory of individual SRA runs as seperate subdirs.
	seq="/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/workdir2/IDFT_NS.fasta"
	prodigal -i $seq -g $i -o ./verb_t_"$i".gff -f "gff" -a ./trans_file_t_"$i".faa -d ./nuc_file_t_"$i".fna
done
