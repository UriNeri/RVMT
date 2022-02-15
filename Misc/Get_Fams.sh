#!/bin/bash
#hostname

Sealsoader
cd /media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/New_ent/

COGS=(COG0741 COG4824 COG4764 COG4678 COG2951 COG3772 COG0860 COG0739 COG4193) # COG3179
mkdir COGs Pfams TIGR IPR CATH FunFams Cdd
for i in $COGS[*]; do
	wget -P ./COGs/ https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fasta/"$i".fa.gz
	wget -P ./COGs/ https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fasta/"$i".tsv.gz
	extract ./COGs/"$i".fa.gz
	extract ./COGs/"$i".tsv.gz
	Rscript /media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/Super_Misc/Trim2Core.r /media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/New_ent2/COGs/"$i".fa /media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/New_ent2/COGs/"$i".tsv

done

Pfam=(PF01520 PF13529 PF01183 PF01476 PF05257 PF01551 PF10464 PF05497 PF05838 PF06737 PF03245 PF05106 PF16754 PF00959 PF10746 PF13702 PF01464 PF18341 PF16079 PF11860 PF04688 PF16082 PF07332 PF11031 PF00144)
for i in $Pfam[*]; do
	# srcl="http://pfam.xfam.org/family/"$i"/alignment/long/gzipped"
	srcl="http://pfam.xfam.org/family/"$i"/alignment/seed/format?format=fasta&alnType=seed&order=t&case=l&gaps=dashes&download=1"
	wget -O ./Pfams/"$i".afa $srcl
	fa_strict ./Pfams/"$i".afa -l=1 >./Pfams/"$i".faa
	# extract ./Pfams/"$i".fa.gz
done

TIGRs=(TIGR03495 TIGR01593 TIGR01606 TIGR01594 TIGR02283 TIGR02282 TIGR02883)
for i in $TIGRs[*]; do
	srcl="https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.SEED/"$i".1.SEED"
	wget -O ./TIGR/"$i".sto $srcl
	reformat.pl sto fas ./TIGR/"$i".sto ./TIGR/"$i".faa
	fa_strict ./TIGR/"$i".faa -l=1 >./TIGR/"$i".fa
done

# IPRs=(IPR031922 IPR033907 IPR001165 IPR002196 IPR034690) #IPR023346

CATHs=(1.10.530.40 1.10.3690.10 1.20.141.10 1.10.10.1120 3.40.630.40)
for i in $CATHs[*]; do
	srcl="ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/sequence-data/sequence-by-superfamily/cath-superfamily-seqs-"$i".fa"
	wget -O ./CATH/"$i".fa $srcl
done

# FunFams in CATHs
CATHs=(1.10.530.40 1.10.3690.10 1.20.141.10 1.10.10.1120 3.40.630.40)
for ix in $CATHs[*]; do
	for ((iy = 0; iy < 100; iy++)); do
		srcl="http://www.cathdb.info/version/v4_3_0/superfamily/"$ix"/funfam/"$iy"/files/stockholm?task_id=&max_sequences=2000&onlyseq=1"
		wget -O ./FunFams2/"$ix".FF."$iy".sto $srcl
	done
done

# FunFams in 1.10.101.10
FunFams=(3 17 10 20 24 26 31 30 36 21 16 14 12)
for i in $FunFams[*]; do
	srcl="http://www.cathdb.info/version/v4_3_0/superfamily/1.10.101.10/funfam/"$i"/files/stockholm?task_id=&max_sequences=20000"
	wget -O ./FunFams/1.10.101.10.FF."$i".sto $srcl
	reformat.pl sto fas ./FunFams/1.10.101.10.FF."$i".sto ./FunFams/1.10.101.10.FF."$i".afa
done

# FunFams in 2.70.70.10
FunFams=(2 3 4 5 6 9 10 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 29 30 31 32 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 55)
for i in $FunFams[*]; do
	srcl="http://www.cathdb.info/version/v4_3_0/superfamily/2.70.70.10/funfam/"$i"/files/stockholm?task_id=&max_sequences=20000"
	wget -O ./FunFams/2.70.70.10.FF."$i".sto $srcl
	reformat.pl sto fas ./FunFams/2.70.70.10.FF."$i".sto ./FunFams/2.70.70.10.FF."$i".afa
done

CDDs=(cd00118 cd16904 cd13402 cd16903 cd16902 cd16901 cd16900 cd16899 cd16889 cd13401 cd12799 cd00737 cd00736 cd00735 cd00442 cd00325 cd00254 cd00119 cd16897 cd02696 cd16888) #cd00978
for i in $CDDs[*]; do
	# srcl="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="$i
	srcl="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="$i"&seqout=1&maxaln=-1"
	wget -O ./Cdd/"$i".afa $srcl
	sed -i '/^</d' ./Cdd/"$i".afa
	fa_strict ./Cdd/"$i".afa -l=1 >./Cdd/"$i".faa
	print $i $srcl
done

###### MISC ######
for i in ./*".FASTA"; do
	ini_name=$(basename $i ".FASTA")
	mv $i "$ini_name"".faa"
done
cat ./*.faa | fa_strict -l=1 >./MegaMerged.fa

Sealsoader
alias Profiler.sh='/media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/Super_Misc/Profiler_motifs.sh '
THREADS=10 #t
output_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/Custom2/"
params='nanabanana' #s
input_fasta="/media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/Custom/MegaMerged.faa"
Memory=9000            #M
min_prec_id=0.99       #d
min_prec_cov=0.95      #c
MCL_inflation=1.8      #f
min_nseq=2             #N
Precluster="True"      #p
Max_sensitivity="True" #x
Cls_Prefix="Cls.Lys"   #P
Use_singlt="True"      #L
SkipCA="False"         #k
AddCons=True           #a
rm_tmp="False"         #r
Use_Mafft="True"       #S
addss="True"
Profiler.sh -t $THREADS -M $Memory -o $output_dir -i $input_fasta -s $params -d $min_prec_id -c $min_prec_cov -f $MCL_inflation -N $min_nseq -p $Precluster -L $Use_singlt -x $Max_sensitivity -P $Cls_Prefix -k $SkipCA -a $AddCons -r $rm_tmp -S $Use_Mafft -H $addss

cat ./*.faa | fa_strict -l=1 >./MegaMerged.fa

for i in ./*".afa"; do
	ini_name=$(basename $i ".afa")
	mv $i "$ini_name".fas
	reformat.pl "$ini_name".fas "$ini_name".sto
done

for i in ./Trimmed_*".fa"; do
	ini_name=$(basename $i ".fa")
	Nini_name=$(echo $ini_name | sed 's|Trimmed_||g')
	# echo $Nini_name
	mafft --thread "$THREADS" --auto $i >"$Nini_name".fas
	reformat.pl "$Nini_name".fas "$Nini_name".sto
done

####################

for i in ./*; do
	ini_name=$(basename $i ".fa")
	ini_name=$(basename $ini_name ".faa")
	ini_name=$(basename $ini_name ".sto")
	ini_name=$(basename $ini_name ".afa")
	ini_name=$(basename $ini_name ".fas")
	seqkit rmdup -s $i -i -o tmp
	mv tmp $i

	mafft --thread "$THREADS" --auto $i >../MSA/"$ini_name".fas

done

cd ../MSA
for i in ./*".fas"; do
	ini_name=$(basename $i ".fas")
	reformat.pl $i "$ini_name".sto
done
#
mv *.sto ../STO/

cd ../STO/
for i in ./*".sto"; do
	ini_name=$(basename $i ".sto")
	addss.pl -i $i -o "$ini_name".a3m -clu
done

mv *.a3m ../A3M/
cd ../A3M/
for i in ./*".a3m"; do
	ini_name=$(basename $i ".a3m")
	hhmake -v 2 -name "$ini_name" -i $i -o ../HHM/"$ini_name".hhm -id 100 -diff -M
done

cd ../MSA/
for i in ./*".fas"; do
	ini_name=$(basename $i ".fas")
	hmmbuild -n "$ini_name" ../HMM/"$ini_name".hmm $i
done

# for i in ./*".sto"; do
#   ini_name=$(basename $i ".sto")
#   reformat.pl  "$ini_name".sto ../Merge_fas/"$ini_name".fas
# done
##2

ffindex_build -s ../hhdb/db_msa.ff{data,index} /media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/New_ent/Merge_fas/
OMP_NUM_THREADS=1
mpirun -np 1 ffindex_apply db_msa.ff{data,index} \
	-i db_a3m_wo_ss.ffindex -d db_a3m_wo_ss.ffdata \
	-- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0
rm db_msa.ff{data,index}

mpirun -np 1 ffindex_apply db_a3m_wo_ss.ff{data,index} \
	-i db_a3m.ffindex -d db_a3m.ffdata -- addss.pl -v 0 stdin stdout
rm db_a3m_wo_ss.ff{data,index}

mpirun -np 1 ffindex_apply db_a3m.ff{data,index} \
	-i db_hhm.ffindex -d db_hhm.ffdata -- hhmake -i stdin -o stdout -v 0

mpirun -np 1 cstranslate -f -x 0.3 -c 4 -I a3m -i db_a3m -o db_cs219

sort -k3 -n -r db_cs219.ffindex | cut -f1 >sorting.dat

ffindex_order sorting.dat db_hhm.ff{data,index} db_hhm_ordered.ff{data,index}
mv db_hhm_ordered.ffindex db_hhm.ffindex
mv db_hhm_ordered.ffdata db_hhm.ffdata

ffindex_order sorting.dat db_a3m.ff{data,index} db_a3m_ordered.ff{data,index}
mv db_a3m_ordered.ffindex db_a3m.ffindex
mv db_a3m_ordered.ffdata db_a3m.ffdata

####### Global enviroment #######
THREADS=8      # uncomment if not passed via the qsub (-v) command
MEMORY=50g     # uncomment if not passed via the qsub (-v) command
algo="hhblits" # uncomment if not passed via the qsub (-v) command (hhsuite sepretaly)
eval=0.00001
scov=1
qcov=50

input=/media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/Custom3/ALL_nuc_3007_sixframe.faa
# input=/media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/Custom3/testout.faa
input_dir=$(dirname $input)
OUTDIR=/media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/Custom3/
profileDir=/media/HDD1/uri/RNA_Vir_MTs/V3/Darius/lys/Custom3/
SDB_path=$profileDir/hhdb/db
ini_name=$(basename $input ".faa")
cd $OUTDIR
mkdir HHspilitted_output spilitted
seqkit split -i $input -O spilitted -j "$THREADS"
cd spilitted
for i in ./*".fas"; do
	ini_name=$(basename $i ".fas")
	# mv $i "$ini_name".fas
	addss.pl -i $i -o ../spilitted_a3m/"$ini_name".a3m -fas
done

ls spilitted >tmparg_file.txt
sed -i 's| |\n|g' tmparg_file.txt
parallel -a tmparg_file.txt -j"$THREADS" hhblits -v 1 -pre_evalue_thresh 1 -maxfilt 200 -n 1 -M 50 -i spilitted/{} -id 100 -cov $qcov -d "$SDB_path" -e "$eval" -o /dev/null -cpu 1 -hide_cons -Z 100000 -blasttab HHspilitted_output/HHMsearch_$(basename {.} .faa)_hhsearch_rawout.tsv
# parallel  -a tmparg_file.txt -j"$THREADS" hhblits -n 2 -M 50 -i spilitted/{} -d "$SDB_path"  -e "$eval" -o HHspilitted_output/HHMsearch_out_$(basename {.} .faa)_hhsearch_rawout.tsv -cpu 1 -hide_cons -Z 100000 -blasttab HHspilitted_output/HHMsearch_$(basename {.} .faa)_hhsearch_rawout.tsv

parallel -a tmparg_file.txt -j"$THREADS" hhsearch -n 1 -M 50 -i spilitted/{} -d "$SDB_path" -e "$eval" -o /dev/null -cpu 1 -hide_cons -Z 100000 -blasttab HHspilitted_output/HHMsearch_$(basename {.} .faa)_hhsearch_rawout.tsv
# parallel  -a tmparg_file.txt -j"$THREADS" hhsearch -n 2 -M 50 -i spilitted/{} -d "$SDB_path"  -e "$eval" -o HHspilitted_output/HHMsearch_out_$(basename {.} .faa)_hhsearch_rawout.tsv -cpu 1 -hide_cons -Z 100000 -blasttab HHspilitted_output/HHMsearch_$(basename {.} .faa)_hhsearch_rawout.tsv
# parallel  -a tmparg_file.txt -j"$THREADS" hhsearch -maxfilt 1000 -pre_evalue_thresh 0.01 -n 1 -M 50 -i "$ini_name"_splitted/{} -d $SDB_path  -e $eval -o /dev/null -cpu 1 -hide_cons -Z 100000 -blasttab HHMsearch_"$ini_name"/$(basename {.} .faa)_hhsearch_rawout.tsv
# parallel  -a tmparg_file.txt -j"$THREADS" hhblits -v 0 -maxfilt 200 -norealign -pre_evalue_thresh 0.001 -n 1 -M 50 -i "$ini_name"_splitted/{} -d $SDB_path  -e $eval -o /dev/null -cpu 1 -hide_cons -Z 100000 -blasttab HHMsearch_"$ini_name"/$(basename {.} .faa)_hhsearch_rawout.tsv
# parallel  -a tmparg_file.txt -j"$THREADS" hhblits -noprefilt -norealign -n 1 -M 50 -i "$ini_name"_splitted/{} -d $SDB_path  -e $eval -o /dev/null -cpu 1 -hide_cons -Z 100000 -blasttab HHMsearch_"$ini_name"/$(basename {.} .faa)_hhsearch_rawout.tsv

#TO-DO: add file name to each line of each file
cat HHMsearch_"$ini_name"/* >"$out_pr"_HHMsearch_vs_"$ini_name"_vs_"$out_su"_rawout.tsv

mkdir /mnt/tmp_hhdb/

mount -t tmpfs -o size=10g tmpfs /mnt/tmp_hhdb/
cp $profileDir/hhdb/ /mnt/tmp_hhdb/hhdb/
SDB_path=/mnt/tmp_hhdb/hhdb/db

mount -t tmpfs -o size=512m tmpfs /mnt/ramdisk

################ MISC ################
for i in ali.210416f/*".fas"; do
	ini_name=$(basename $i ".fas")
	reformat.pl fas a2m "$ini_name".fas -o ../a3mfiles_woss/"$ini_name".a2m
done

for i in ./*".fas"; do
	ini_name=$(basename $i ".fas")
	reformat.pl fas a2m "$ini_name".fas -o ../a3mfiles_woss/"$ini_name".a2m
done

for i in ./*".sto"; do
	ini_name=$(basename $i ".sto")
	addss.pl -i $i -o ../a3mfiles/"$ini_name".a3m -clu
done

for i in ./*".a3m"; do
	ini_name=$(basename $i ".a3m")
	hhmake -v 2 -name "$ini_name" -i $i -o ../Merge_hhm/"$ini_name".hhm -id 100 -diff -M
done

for i in ./*".a2m"; do
	ini_name=$(basename $i ".a2m")
	hhmake -v 2 -name "$ini_name" -i $i -o ../hhfiles/"$ini_name".hhm -id 100 -diff -M
done

for i in ./*".sto"; do
	ini_name=$(basename $i ".sto")
	hhmake -v 2 -name "$ini_name" -i $i -o ../hhfiles/"$ini_name".hhm -id 100 -diff -M
done

for i in ./*".hhm"; do
	ini_name=$(basename $i ".hhm")
	addss.pl $i ../hhfiles_wss/"$i" -hmm
done
