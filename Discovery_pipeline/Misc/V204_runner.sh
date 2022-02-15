module load R/R-3.3.2
THREADS=48
Neri_scripts=/urigo/urineri/scripts/RNA_Viruses/Neri_scripts/
resources_dir="$Neri_scripts"/resources
V2_dir=/scratch200/urineri/V2/

step=4
workdir="$V2_dir"/MTs/step_"$step"/
refdirs="$V2_dir"/MGs/
cd $workdir
Rscript /urigo/urineri/scripts/RNA_Viruses/Neri_scripts/presence_absence/Rscripts/V204_presence_absence_sys_args.r
