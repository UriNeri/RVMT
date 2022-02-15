#!/bin/bash
#hostname
###### Submit a "Iter_RdRp_Search_psihmmseqs2" jobs  ######
cd /urigo/urineri/logs/
Set_ID=0

rpD=40 # RAM in GB per Diamond BLASTp
tpD=10 #Threads per Diamond BLASTp

rpM=350 # RAM in GB per mmseqs2
tpM=46  #Threads per mmseqs2

rpP=58 # RAM in GB per PSI-BALST
tpP=12 #Threads per PSI-BALST

rpH=350 # RAM in GB per HMMsearch
tpH=46  #Threads per HMMsearch

qsub -v set_ID="$Set_ID",THREADS="$tpD",algo="diamond" -l mem="$rpD"gb,nodes=1:ppn=$tpD -q hugemem /urigo/urineri/scripts/RNA_Viruses/Neri_scripts/RdRp_iterative_search/V204_runner_Iter_RdRp_Search_psihmmseqs2.sh -m be -M urineri@mail.tau.ac.il
qsub -v set_ID="$Set_ID",THREADS="$tpH",algo="hmmsearch" -l mem="$rpH"gb,nodes=1:ppn=$tpH -q gophnabig /urigo/urineri/scripts/RNA_Viruses/Neri_scripts/RdRp_iterative_search/V204_runner_Iter_RdRp_Search_psihmmseqs2.sh -m be -M urineri@mail.tau.ac.il
qsub -v set_ID="$Set_ID",THREADS="$tpM",algo="mmseqs2" -l mem="$rpM"gb,nodes=1:ppn=$tpM -q gophnabig /urigo/urineri/scripts/RNA_Viruses/Neri_scripts/RdRp_iterative_search/V204_runner_Iter_RdRp_Search_psihmmseqs2.sh -m be -M urineri@mail.tau.ac.il
qsub -v set_ID="$Set_ID",THREADS="$tpP",algo="psiblast" -l mem="$rpP"gb,nodes=1:ppn=$tpP -q hugemem /urigo/urineri/scripts/RNA_Viruses/Neri_scripts/RdRp_iterative_search/V204_runner_Iter_RdRp_Search_psihmmseqs2.sh -m be -M urineri@mail.tau.ac.il

###### Submit a "HHtools4CS.sh" job  ######

rpH=350 # RAM in GB per HMMsearch
tpH=46  #Threads per HMMsearch
qsub -v set_ID="$Set_ID",THREADS="$tpH",Memory="$rpH",db_name="pfam_32" -l mem="$rpH"gb,nodes=1:ppn=$tpH -q gophnabig /urigo/urineri/scripts/RNA_Viruses/Neri_scripts/RdRp_iterative_search/HHtools4CS.sh -m be -M urineri@mail.tau.ac.il
