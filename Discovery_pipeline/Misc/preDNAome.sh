#!/bin/bash
#hostname

# module load MMseqs2/Jan2020-avx2 # == MMseqs2 7th Jan commit
# export PATH=$PATH:/urigo/urineri/resources/mmseqs/bin/

# Neri_scripts=/media/uri/HDD1/uri/V1_hard_mode/V2/Neri_scripts/
Neri_scripts=/urigo/urineri/scripts/RNA_Viruses/Neri_scripts/
resources_dir="$Neri_scripts"/resources
mmseqs2="$resources_dir"/mmseqs_SSE4/bin/mmseqs # mmseqs2_SSE4
# mmseqs2="$resources_dir"/mmseqs_AVX2/bin/mmseqs # mmseqs2_AVX2
# export PATH=$PATH:"$resources_dir"/mmseqs_SSE4/bin/
export PATH=$PATH:"$resources_dir"/mmseqs_AVX2/bin/
# export PATH=$PATH:/urigo/DBs/resources/mmseqs/bin/
THREADS=48

V2_dir=/scratch200/urineri/V2/
V2_dir=/davidb/urineri/V2/

MGs="$V2_dir""MGs/"
cd $MGs
# cd ./s2/
# i="s2"

# study_list=("Grasslands_soil_microbial_communities_from_the_Angelo_Coastal_Reserve__California__USA_1")  #  "Enriched_cultures_of_PCE_dechlorinating_microbial_communities_from_Ithaca__New_York__USA" "Enriched_soil_aggregate_microbial_communities_from_Iowa_State_university_to_study_microbial_drivers_of_carbon_cycling" "Freshwater_microbial_communities_amended_with_dissolved_organic_matter__DOM__from_various_rivers_in_the_United_States" "Freshwater_sediment_methanotrophic_microbial_communities_from_Lake_Washington_under_simulated_oxygen_tension" "Lab_enriched_sorghum-adapted_microbial_communities_from_Joint_BioEnergy_Institute__Emeryville__California__United_States" "Marine_archaeal_communities_from_Monterey_Bay__CA__that_are_ammonia_oxidizing" "Marine_microbial_communities_from_expanding_oxygen_minimum_zones_in_the_northeastern_subarctic_Pacific_Ocean" "Marine_microbial_communities_from_the_Southern_Atlantic_ocean_transect_to_study_dissolved_organic_matter_and_carbon_cycling_MTs_only_low_euk")
# for i in ${study_list[*]};
# do

# study_name=$(basename $i)
# mkdir "$i"/MMseqs2DB
# cd "$i"/MMseqs2DB/
# echo "Started easy-linclust for $study_name"
# mmseqs easy-linclust ../sk100.fna sk100.clu tmp --min-seq-id 0.99 -c 1.0 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS --compressed 1
# mv ./sk100.clu_rep_seq.fasta ../sk100_reps.fna
# cd ../
# rm ./MMseqs2DB/ -r

# echo "Finished easy-linclust for $study_name"
# cd ../

# done
# cd ../

# for j in ("s3" "s3")
# do
# cd "$MGs"/"$j"/

# for i in./*
# do

# study_name=$(basename $i)
# mkdir "$i"/MMseqs2DB
# cd "$i"/MMseqs2DB/
# echo "Started easy-linclust for $study_name"
# mmseqs easy-linclust ../sk100.fna sk100.clu tmp --min-seq-id 0.99 -c 1.0 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS --compressed 1
# mv ./sk100.clu_rep_seq.fasta ../sk100_reps.fna
# cd ../
# rm ./MMseqs2DB/ -r

# echo "Finished easy-linclust for $study_name"
# cd ../

# done
# echo "Finished easy-linclust for $j"

# done

# db_path="$V2_dir"IMG_VR/MMseqs2DB/
# cd $db_path
# mmseqs easy-linclust ../linear_IMGVR_DNA_only_scafs.fna DNA_virDB.clu tmp --min-seq-id 0.99 -c 1.0 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS

# for i in ./*/
# done
# namma=$(basename $i)
# cp2="/run/user/1000/gvfs/sftp:host=powerlogin-be1.tau.ac.il,user=urineri/scratch200/urineri/V2/MGs/s1/""$namma"
# cp "$i""/sk100.fna" "$cp2"/sk100.fna
# done

# ##### For prots #####
# # THREADS=8
# mkdir MMseqs2DB
# cd MMseqs2DB
# mmseqs createdb ../dmnd.faa DB
# mmseqs linclust  DB DB_clu tmp --min-seq-id 0.95 -c 1 --cov-mode 1 --kmer-per-seq-scale 0.5  --threads $THREADS
# mmseqs result2repseq DB DB_clu DB_clu_rep --threads $THREADS
# mmseqs result2flat DB DB DB_clu_rep DB_clu_rep.fasta  --use-fasta-header
# mv DB_clu_rep.fasta ../dmnd_reps.faa

# mkdir MMseqs2DB
# cd MMseqs2DB
# mmseqs createdb ../cated_prots.faa  DB
# mmseqs linclust  DB DB_clu tmp --min-seq-id 0.95 -c 1 --cov-mode 1 --kmer-per-seq-scale 0.5  --threads $THREADS
# mmseqs result2repseq DB DB_clu DB_clu_rep --threads $THREADS
# mmseqs result2flat DB DB DB_clu_rep DB_clu_rep.fasta  --use-fasta-header
# mv DB_clu_rep.fasta ../cated_prots_reps.faa
# cd ../
# rm MMseqs2DB -r

##### NEW #####
# study_list=("s1" "s2")
# for i in ${study_list[*]};
# do
# cd $MGs/"$i"
# mkdir sk100
# cd  sk100
# echo "Started easy-linclust for $i"
# mmseqs easy-linclust ../cated_sk100_reps.fna sk100.clu tmp --min-seq-id 0.99 -c 1.0 --cov-mode 1 --kmer-per-seq-scale 0.4  --threads $THREADS --compressed 1
# # mv ./sk100.clu_rep_seq.fasta ../sk100_reps.fna
# echo "Finished easy-linclust for $study_name"
# cd $MGs
# done

# study_list=("s1" "s2" "s3") # "s2")
study_list=("s2" "s3") # "s2")
study_list=("s2")

for i in ${study_list[*]}; do
	cd $MGs/"$i"
	mkdir MMseqs2DB
	echo "Started formating MMseqs2DB for $i"
	mmseqs createdb ./cated_sk100_reps.fna ./MMseqs2DB/MMseqs2DB --dbtype 0 # --compressed 1
	# mmseqs  createindex ./MMseqs2DB/MMseqs2DB ./tmp  --threads $THREADS  --search-type 3
	echo "Finished formating MMseqs2DB for $i"
	# cd $MGs
done

# cp /scratch200/urineri/V2/MGs/s3/head_cated_sk100_reps.fna /a/home/cc/students/lifesci/urineri/public_html/Data/4Milot/head_cated_sk100_reps.fna
# awk '{if(NR%2==0)print}' head_cated_sk100_reps.fna > awked_head_cated_sk100_reps.seqs

# mmseqs createdb ./head_cated_sk100_reps.fna ./MMseqs2DB/MMseqs2DB    --dbtype 0 #--compressed 1
# mmseqs  createindex ./MMseqs2DB/MMseqs2DB ./tmp  --threads $THREADS  --search-type 3 --compressed 1
