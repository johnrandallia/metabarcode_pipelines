#! /bin/bash

# Script to run metamate

# Load conda environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate metamate_env

######################################################################################################################################################################################
# 
# pre denoising
# 
#######################################################################################################################################################################################

# no denoiser minAbund 2
metamate find \
 -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/UNIQUESbyLIB/all_reads_UNIQUESbyLIB_minAbund2_onlyArthropods.fasta \
 -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/UNIQUESbyLIB/ASV_table_minAbund2_onlyArthro.csv \
 -S /Datos/tools/metamate/specifications_total_and_lib_extensive_tax_binning_no_clade.txt \
 -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta --expectedlength 418  \
 --percentvar 0 --table 5 -o metamate_without_denoiser_minAbund2 -t 18 --overwrite


######################################################################################################################################################################################
# 
# UNOISE3
# 
#######################################################################################################################################################################################
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/UNOISE_relabelled*/unoise3_uchime_onlyArthropods.fasta) ; do
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "Metamate is executed for unoise with alpha $alpha .."
    
    if [ $alpha == "a4" ] || [ $alpha == "a3" ] || [ $alpha == "a2" ] || [ $alpha == "a1" ]
    then
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive.txt"
    else
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive_tax_binning_no_clade.txt"
    fi

    if [ $parentDir != "UNOISE_relabelled_joined" ]
    then
        if [[ $parentDir != *"minsize"* ]]
        then
            metamate find -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/unoise3_uchime_onlyArthropods.fasta \
            -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/ASV_table_unoise3.csv \
            -S $specFile -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta \
            --expectedlength 418  --percentvar 0 --table 5 -o metamate_unoise_$alpha -t 18 --overwrite
            
        fi
    fi
done

# unoise + lulu
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/UNOISE_relabelled*/unoise3_lulu.fasta) ; do
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "Metamate is executed for unoise with alpha $alpha .."
    
    if [ $alpha == "a4" ] || [ $alpha == "a3" ] || [ $alpha == "a2" ] || [ $alpha == "a1" ]
    then
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive.txt"
    else
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive_tax_binning_no_clade.txt"
    fi

    if [ $parentDir != "UNOISE_relabelled_joined" ]
    then
        if [[ $parentDir != *"minsize"* ]]
        then
            metamate find -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/unoise3_lulu.fasta \
            -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/lulu_curated_table.txt \
            -S $specFile -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta \
            --expectedlength 418  --percentvar 0 --table 5 -o metamate_unoise_lulu_$alpha \
            -t 18 --overwrite 
        fi
    fi
done


######################################################################################################################################################################################
# 
# DADA 2
# 
#######################################################################################################################################################################################
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2_omega*/dada2_onlyArthropods.fasta) ; do
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    omega="$(sed -E 's/dada2_omega//g' <(echo $parentDir))"
    # echo $omega
    echo "Metamate is executed for dada2 with omega_c $omega .."

    specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive.txt"

    metamate find -A $D \
    -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/ASV_table_dada2.csv \
    -S $specFile -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta \
    --expectedlength 418  --percentvar 0 --table 5 -o metamate_dada2_$omega -t 18 --overwrite --realign
done

for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2_omega*/dada2_lulu.fasta) ; do
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    omega="$(sed -E 's/dada2_omega//g' <(echo $parentDir))"
    # echo $omega
    echo "Metamate is executed for dada2 with omega_c $omega .."

    specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive.txt"

    metamate find -A $D \
    -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/lulu_curated_table.txt \
    -S $specFile -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta \
    --expectedlength 418  --percentvar 0 --table 5 -o metamate_dada2_lulu_$omega -t 18 --overwrite --realign
done


######################################################################################################################################################################################
# 
# DnoisE
# 
#######################################################################################################################################################################################
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dnoise_relabelled*/dnoise_uchime_onlyArthropods.fasta) ; do
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "Metamate is executed for dnoise with alpha $alpha .."
    
    if [ $alpha == "alpha1" ]
    then
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive.txt"
    else
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive_tax_binning_no_clade.txt"
    fi

    metamate find -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/dnoise_uchime_onlyArthropods.fasta \
    -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/ASV_table_dnoise.csv \
    -S $specFile -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta \
    --expectedlength 418  --percentvar 0 --table 5 -o metamate_dnoise_$alpha -t 18 --overwrite 

done

for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dnoise_relabelled*/dnoise_lulu.fasta) ; do
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "Metamate is executed for dnoise with alpha $alpha .."
    
    if [ $alpha == "alpha1" ]
    then
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive.txt"
    else
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive_tax_binning_no_clade.txt"
    fi

    metamate find -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/dnoise_lulu.fasta \
    -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/lulu_curated_table.txt \
    -S $specFile -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta \
    --expectedlength 418  --percentvar 0 --table 5 -o metamate_dnoise_lulu_$alpha -t 18 --overwrite --realign

done


######################################################################################################################################################################################
# 
# SWARM
# 
#######################################################################################################################################################################################
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d*/swarm_onlyArthropods.fasta) ; do
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "Metamate is executed for swarm with d $alpha .."
    
    if [ $alpha == "d10" ]
    then
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive.txt"
    else
        specFile="/Datos/tools/metamate/specifications_total_and_lib_extensive_tax_binning_no_clade.txt"
    fi

    metamate find -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/swarm_onlyArthropods.fasta \
    -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/ASV_table_swarm.csv \
    -S $specFile -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta \
    --expectedlength 418  --percentvar 0 --table 5 -o metamate_swarm_$alpha -t 18 --overwrite --realign

done

# swarm lulu
for D in $(find /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d*/swarm_lulu_onlyArthropods.fasta) ; do
    parentPath="$(dirname "$D")"
    parentDir="$(basename $parentPath)"
    # echo $parentDir
    alpha="$(sed -E 's/.*_//g' <(echo $parentDir))"
    # echo $alpha
    echo "Metamate is executed for swarm+lulu with d $alpha .."
    
    metamate find -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/swarm_lulu_onlyArthropods.fasta \
    -M /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/$parentDir/lulu_curated_table.txt \
    -S /Datos/tools/metamate/specifications_total_and_lib_extensive_tax_binning_no_clade.txt -R /Datos/mockcommunities_EPA45/EPA45_chosen_418bp_derep_stopremoved_geneious_uniqHeader_minLen410_onlyArthropods.fasta \
    --expectedlength 418  --percentvar 0 --table 5 -o metamate_swarm_lulu_$alpha -t 18 --overwrite --realign

done



#################################################################################################################################
#
# METAMATE DUMP
#
#################################################################################################################################

# < 5p vnaASV
# no denoiser
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/UNIQUESbyLIB/all_reads_UNIQUESbyLIB_minAbund2_onlyArthropods.fasta -C ../resultcache -o metamate_no_denoiser_filtered --overwrite -i 328
mafft --thread 15 --retree 1 metamate_no_denoiser_filtered.fasta > metamate_no_denoiser_filtered_aligned.fasta

# unoise3
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/UNOISE_relabelled_joined_a2/unoise3_uchime_onlyArthropods.fasta -C ../resultcache -o metamate_unoise3_filtered --overwrite -i 1251
mafft --thread 15 --retree 1 metamate_unoise3_filtered.fasta > metamate_unoise3_filtered_aligned.fasta

# unoise + lulu
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/UNOISE_relabelled_joined_a2/unoise3_lulu.fasta -C ../resultcache -o metamate_unoise3_lulu_filtered --overwrite -i 1455
mafft --thread 15 --retree 1 metamate_unoise3_lulu_filtered.fasta > metamate_unoise3_lulu_filtered_aligned.fasta

# dnoise
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dnoise_relabelled_joined_alpha5/dnoise_uchime_onlyArthropods.fasta -C ../resultcache -o metamate_dnoise_filtered --overwrite -i 328
mafft --thread 15 --retree 1 metamate_dnoise_filtered.fasta > metamate_dnoise_filtered_aligned.fasta

# dnoise + lulu
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dnoise_relabelled_joined_alpha5/dnoise_lulu.fasta -C ../resultcache -o metamate_dnoise_lulu_filtered --overwrite -i 1379
mafft --thread 15 --retree 1 metamate_dnoise_lulu_filtered.fasta > metamate_dnoise_lulu_filtered_aligned.fasta

# dada2
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2/dada2_onlyArthropods.fasta -C ../resultcache -o metamate_dada2_filtered --overwrite -i 1354
mafft --thread 15 --retree 1 metamate_dada2_filtered.fasta > metamate_dada2_filtered_aligned.fasta

# dada2 + lulu
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2/dada2_lulu.fasta -C ../resultcache -o metamate_dada2_lulu_filtered --overwrite -i 1456
mafft --thread 15 --retree 1 metamate_dada2_lulu_filtered.fasta > metamate_dada2_lulu_filtered_aligned.fasta

# swarm 
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d1/swarm_onlyArthropods.fasta -C ../resultcache -o metamate_swarm_filtered --overwrite -i 3
mafft --thread 15 --retree 1 metamate_swarm_filtered.fasta > metamate_swarm_filtered_aligned.fasta

# swarm + lulu
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d1/swarm_lulu_onlyArthropods.fasta -C ../resultcache -o metamate_swarm_lulu_filtered --overwrite -i 337
mafft --thread 15 --retree 1 metamate_swarm_lulu_filtered.fasta > metamate_swarm_lulu_filtered_aligned.fasta

python3 /Datos/tools/metamate/metamate/metamate.py dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d1/swarm_lulu_onlyArthropods.fasta -C ../resultcache -o metamate_swarm_lulu_filtered --overwrite -i 337




# < 5p estimated naASVs
# no denoiser
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/UNIQUESbyLIB/all_reads_UNIQUESbyLIB_minAbund2_onlyArthropods.fasta -C ../resultcache -o metamate_no_denoiser_filtered --overwrite -i 387
mafft --thread 15 --retree 1 metamate_no_denoiser_filtered.fasta > metamate_no_denoiser_filtered_aligned.fasta

# unoise3
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/UNOISE_relabelled_joined_a2/unoise3_uchime_onlyArthropods.fasta -C ../resultcache -o metamate_unoise3_filtered --overwrite -i 1359
mafft --thread 15 --retree 1 metamate_unoise3_filtered.fasta > metamate_unoise3_filtered_aligned.fasta

# unoise + lulu
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/UNOISE_relabelled_joined_a2/unoise3_lulu.fasta -C ../resultcache -o metamate_unoise3_lulu_filtered --overwrite -i 1555
mafft --thread 15 --retree 1 metamate_unoise3_lulu_filtered.fasta > metamate_unoise3_lulu_filtered_aligned.fasta

# dnoise
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dnoise_relabelled_joined_alpha5/dnoise_uchime_onlyArthropods.fasta -C ../resultcache -o metamate_dnoise_filtered --overwrite -i 1626
mafft --thread 15 --retree 1 metamate_dnoise_filtered.fasta > metamate_dnoise_filtered_aligned.fasta

# dnoise + lulu
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dnoise_relabelled_joined_alpha5/dnoise_lulu.fasta -C ../resultcache -o metamate_dnoise_lulu_filtered --overwrite -i 358
mafft --thread 15 --retree 1 metamate_dnoise_lulu_filtered.fasta > metamate_dnoise_lulu_filtered_aligned.fasta

# dada2
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2/dada2_onlyArthropods.fasta -C ../resultcache -o metamate_dada2_filtered --overwrite -i 1366
mafft --thread 15 --retree 1 metamate_dada2_filtered.fasta > metamate_dada2_filtered_aligned.fasta

# dada2 + lulu
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2/dada2_lulu.fasta -C ../resultcache -o metamate_dada2_lulu_filtered --overwrite -i 1353
mafft --thread 15 --retree 1 metamate_dada2_lulu_filtered.fasta > metamate_dada2_lulu_filtered_aligned.fasta

# swarm 
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d1/swarm_onlyArthropods.fasta -C ../resultcache -o metamate_swarm_filtered --overwrite -i 343
mafft --thread 15 --retree 1 metamate_swarm_filtered.fasta > metamate_swarm_filtered_aligned.fasta

# swarm + lulu
metamate dump -A /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/swarm_relabelled_d1/swarm_lulu_onlyArthropods.fasta -C ../resultcache -o metamate_swarm_lulu_filtered --overwrite -i 378
mafft --thread 15 --retree 1 metamate_swarm_lulu_filtered.fasta > metamate_swarm_lulu_filtered_aligned.fasta

