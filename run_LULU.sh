#! /bin/bash

# script to run LULU
# M. Eisele, C. Andujar, B. Emerson
# July 2024

# Load conda environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate metabarcode_data_analysis_env

################################################################################################################################
#
# UNOISE
#
#################################################################################################################################
N=10
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/UNOISE_relabelled_joined_a*/unoise3_uchime_onlyArthropods.fasta) ;
do
    ((i=i%N)); ((i++==0)) && wait
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "LULU is executed for UNOISE3 with d $alpha .."

    Rscript --vanilla /Datos/tools/LULU.R /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/match_list.txt \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/ASV_table_unoise3.csv \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/unoise3_uchime_onlyArthropods.fasta \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir &

done

for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/UNOISE_relabelled_joined_a*/unoise3_uchime_onlyArthropods.fasta ) ;
do
    echo "running seqkit to extract new fasta from lulu runs for $D.."
    
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"

    seqkit grep -v -n -f /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/lulu_discarded_otus.txt -o /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/unoise3_lulu.fasta $D

done


################################################################################################################################
#
# DADA2
#
#################################################################################################################################
N=10
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2_omega*/dada2_onlyArthropods.fasta) ;
do
    ((i=i%N)); ((i++==0)) && wait
    if [[ ${D} != *"chimeras"* ]];
    then
        parentPath="$(dirname "$D")"
        parentDir="$(basename $parentPath)"
        # echo $parentDir
        alpha="$(sed -E 's/.*_omega//g' <(echo $parentDir))"

        # echo $alpha
        echo "LULU is executed for DADA2 with omega_a $alpha .."  

        Rscript --vanilla /Datos/tools/LULU.R /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/match_list.txt \
        /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/ASV_table_dada2.csv \
        /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/dada2_onlyArthropods.fasta \
        /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir &

    fi
done

for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2_omega*/dada2_onlyArthropods.fasta) ;
do
    if [[ ${D} != *"chimeras"* ]] && [[ ${D} != *"lulu"* ]] ;
    then
        echo "running seqkit to extract new fasta from lulu runs for $D.."
        
        parentPath="$(dirname "$D")"
        parentDir="$(basename $parentPath)"

        seqkit grep -v -n -f /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/lulu_discarded_otus.txt -o /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/dada2_lulu.fasta $D
    fi
done


################################################################################################################################
#
# DnoisE
#
#################################################################################################################################
N=10
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dnoise_relabelled_joined_a*/dnoise_uchime_onlyArthropods.fasta) ;
do
    ((i=i%N)); ((i++==0)) && wait
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "LULU is executed for DnoisE with alpha = $alpha .."

    Rscript --vanilla /Datos/tools/LULU.R /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/match_list.txt \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/ASV_table_dnoise.csv \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/dnoise_uchime_onlyArthropods.fasta \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir &

done

for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dnoise_relabelled_joined_a*/dnoise_uchime_onlyArthropods.fasta) ;
do
    echo "running seqkit to extract new fasta from lulu runs for $D.."
    
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"

    seqkit grep -v -n -f /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/lulu_discarded_otus.txt -o /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/dnoise_lulu.fasta $D

done


################################################################################################################################
#
# SWARM
#
#################################################################################################################################
N=10
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d*/swarm_onlyArthropods.fasta) ;
do
    ((i=i%N)); ((i++==0)) && wait
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "LULU is executed for swarm with d $alpha .."

    Rscript --vanilla /Datos/tools/LULU.R /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/match_list.txt \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/ASV_table_swarm.csv \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/swarm_onlyArthropods.fasta \
    /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir &

done

for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d*/swarm_onlyArthropods.fasta) ; 
do
    echo "running seqkit to extract new fasta from lulu runs.."
    
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"

    seqkit grep -v -n -f /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/lulu_discarded_otus.txt -o /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/swarm_lulu_onlyArthropods.fasta $D

done

