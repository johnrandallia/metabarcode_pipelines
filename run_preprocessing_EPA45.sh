#!/bin/bash

# Script to preprocess raw metabarcoding data

# Load conda environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate metabarcode_data_analysis_env

# Change working directory to the correct location
#cd "$1"
cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed

# Run fastqc on all fastq files
if [ ! -d "fastqc_raw/" ]
then
    echo "Running fastqc on raw reads .."
    mkdir fastqc_raw
    fastqc -o fastqc_raw ../01_raw_Rep1_r1_r2_ordered/*/*.fastq.gz
    echo "Finished running fastqc on raw reads.."
fi


# Trim primers
if [ ! -d "cutadapt/" ]
then
    echo "Running cutadapt..."
    mkdir cutadapt
    cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/01_raw_Rep1_r1_r2_ordered/
    for file in */*.R1.fastq.gz
    do        
        echo ${file}
	    part1="$(echo ${file} | sed 's/\(.*\)R1.*/\1/')" 
	    part2="$(echo ${file} | sed 's/.*R1\(.*\)/\1/')"
	    
	    f2="$part1"
        f2+="R2"
		f2+="$part2"

        basename1="$(basename ${file})"
        r1_output="${basename1%.*}"
        basename2="$(basename ${f2})"
        r2_output="${basename2%.*}"          
 	    
        cutadapt -j 16 -e 0 --untrimmed-output ../02_preprocessed/cutadapt/${r1_output}.no_adapter.fq --untrimmed-paired-output ../02_preprocessed/cutadapt/${r2_output}.no_adapter.fq \
        -g ^CCNGAYATRGCNTTYCCNCG -g ^CNGAYATRGCNTTYCCNCG -g ^NGAYATRGCNTTYCCNCG -g ^GAYATRGCNTTYCCNCG -g ^AYATRGCNTTYCCNCG -G ^TCNGGRTGNCCRAARAAYCA -G ^CNGGRTGNCCRAARAAYCA -G ^NGGRTGNCCRAARAAYCA -G ^GGRTGNCCRAARAAYCA \
        -G ^GRTGNCCRAARAAYCA ${file} ${f2} -o ../02_preprocessed/cutadapt/${r1_output}.cutadapt.fq -p ../02_preprocessed/cutadapt/${r2_output}.cutadapt.fq
        echo "Finished running cutadapt"
    done
fi


# Remove ends with low quality
if [ ! -d "/Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/TRIMMO" ]
then
    echo "Running trimmomatic.."
    cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/cutadapt
    mkdir ../TRIMMO   
    for file in *cutadapt.fq
    do
        trimmomatic SE -threads 16 -phred33 -trimlog ../TRIMMO/$file.log -summary ../TRIMMO/$file.summary.log $file ../TRIMMO/$file.trimmo.fq TRAILING:20
    done
    echo "Finished running trimmomatic.."
fi


# perform another fastqc with trimmed reads
cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/cutadapt
if [ ! -d "../fastqc_after_trimmo" ]
then
    echo "Running fastqc on trimmed reads.."
    cd ../TRIMMO
    mkdir ../fastqc_after_trimmo
    fastqc -o ../fastqc_after_trimmo *.fq
    echo "Finished running fastqc on trimmed reads."
fi


# Pair fastq reads
if [ ! -d "../PAIRED" ]
then
    echo "Pairing reads.."
    cd ../TRIMMO
    mkdir ../PAIRED

    N=16
    for f1 in *R1*trimmo.fq
    do
        ((i=i%N)); ((i++==0)) && wait	    
        echo $f1
	    part1="$(echo $f1 | sed 's/\(.*\)R1.*/\1/')" 
	    part2="$(echo $f1 | sed 's/.*R1\(.*\)/\1/')"
	    
	    f2="$part1"
        f2+="R2"
		f2+="$part2"
		        
	    echo "file1: $f1"
	    echo "file2: $f2"
	    echo "executing pairfq_lite ..."

        fastq_pair $f1 $f2 & 
	    echo "done"
     	
    done
    mv *paired.fq ../PAIRED/
    mv *single.fq ../PAIRED/
fi


# Merge the pairs
if [ ! -d "../MergePair" ]
then
    cd ../PAIRED
    mkdir ../MergePair

    for f1 in *R1*trimmo.fq.paired.fq
    do
	    echo $f1
	    #para inputs
	    part1="$(echo $f1 | sed 's/\(.*\)R1.*/\1/')" 
	    part2="$(echo $f1 | sed 's/.*R1\(.*\)/\1/')"
	    
	    f2="$part1"
	    f2+="R2"
	    f2+="$part2"
	    
        #para outputs
	    part3="$(echo $f1 | sed 's/\(.*\)_Rep1.*/\1/')"
	    #part4="$(echo $f1 | sed 's/\(.*\)_Rep1.R1\(.*\)Ptrim.*/\2/')"
	    part4="_Rep1_"

	    f3="$part3"
	    f3+="$part4"
	    f3+="R1_R2.Ptrim.trimmo.fastq"
	    
	    f4="$part3"
	    f4+="$part4"
	    f4+="R1_R2.Ptrim.trimmo.report.log"
	    
	    echo "file1: $f1"
	    echo "file2: $f2"
	    echo "file3: $f3"

	    echo "executing fastq_mergepairs ..."
	    
	    /Datos/tools/usearch11.0.667_i86linux32 -threads 16 -fastq_mergepairs $f1 -reverse $f2 -fastq_maxdiffs 10 -fastq_pctid 80 -fastqout ../MergePair/$f3 -report ../MergePair/$f4 
	    echo "Merging reads completed."
     done
fi


# discard sequences with a maximum expected error (maxee)>2
if [ ! -d "../maxee1" ]
then
    echo "Running maxee.."  
    cd ../MergePair
    mkdir ../maxee1

    for f1 in *R1_R2.Ptrim.trimmo.fastq
    do
	    echo $f1
	    part1="$(echo $f1 | sed 's/\(.*\)fastq*/\1/')"    

	    f2="$part1"
        f2+="Maxee1.fas"

        f2_fq="$part1"
        f2_fq+="Maxee1.fq"

		f3="$part1"
        f3+="Maxee1.log"

	    echo "input: $f1"
	    echo "output: $f2"
	    echo "logfile: $f3"

	    echo "executing maxee 1 ..."

	    /Datos/tools/usearch11.0.667_i86linux32 -threads 16 -fastq_filter $f1 -fastq_maxee 1 -fastaout ../maxee1/$f2 -fastqout ../maxee1/$f2_fq &> ../maxee1/$f3
	    
	    echo "Finished running maxee."
     done
fi


#LABELLING by library and primer
if [ ! -d "../LABEL_lib" ]
then
    echo "Labelling reads.."
    cd ../maxee1
    mkdir ../LABEL_lib

    for file in *.fas
    do
	    part1="$(echo $file | sed 's/\(.*\)_R1.*.*/\1/')" 
	    label="$part1"
		    
	    echo $file
	    echo $label
        sed "-es/^>\(.*\)/>\1;barcodelabel=$label;/" < ${file} > ../LABEL_lib/${file%.*}.labelbyLIB.fas
    done

    for file in *.fq
    do
	    part1="$(echo $file | sed 's/\(.*\)_R1.*.*/\1/')" 
	    label="$part1"
		    
	    echo $file
	    echo $label
        sed "-es/^@\(.*\)/@\1;barcodelabel=$label;/" < ${file} > ../LABEL_lib/${file%.*}.labelbyLIB.fq
    done
    echo "Finished labelling reads."
fi


# sort by length
# and keeping only reads of 418 (exact length) + - 2 (416-420) on labelling by lib and primer (normal path)
if [ ! -d "../LENGTH416_420" ]
then
    echo "Executing sortbylength 420 with fastas..."
    cd ../LABEL_lib
    mkdir ../LENGTH416_420

    for f1 in *labelbyLIB.fas
    do
	    
	    echo $f1
	    part1="$(echo $f1 | sed 's/\(.*\)fas*/\1/')" 

	    f2="$part1"
        f2+="sorted416_420.fas"
	    f2log="$part1"
	    f2log+="sorted416_420.log"
				    
        echo "input: $f1, output: $f2"
	    
	    /Datos/tools/usearch11.0.667_i86linux32 -threads 16 -sortbylength $f1  -fastaout ../LENGTH416_420/$f2 -minseqlength 416 -maxseqlength 420 &> ../LENGTH416_420/$f2log
	    
	    echo "done sortbylength 416-420 ..."
    done	

    cat ../LENGTH416_420/*.fas >> ../LENGTH416_420/all_reads_LENGTH416_420.fasta

    echo "executing sortbylength 420 with fastqs..."
    for f1 in *labelbyLIB.fq
    do
	    echo $f1
	    part1="$(echo $f1 | sed 's/\(.*\)fq*/\1/')" 

	    f2="$part1"
        f2+="sorted416_420.fq"
				    
        echo "input: $f1, output: $f2"
	    
	    /Datos/tools/usearch11.0.667_i86linux32 -threads 16 -sortbylength $f1  -fastqout ../LENGTH416_420/$f2 -minseqlength 416 -maxseqlength 420
	    
	    echo "done sortbylength 416-420 ..."
    done		

fi


# keep unique (dereplication)
if [ ! -d "../UNIQUESbyLIB" ]
then
    echo "executing fastx_uniques of fastas..."
    cd ../LENGTH416_420
    mkdir ../UNIQUESbyLIB

    for f1 in *sorted416_420.fas
    do
	    
	    #echo $f1
	    part1="$(echo $f1 | sed 's/\(.*\)fas*/\1/')" 

	    f2="$part1"
        f2+="uniques.fas"
        f2_fq="$part1"
	    f2_fq+="uniques.fq"
        f2log="$part1"
	    f2log+="uniques.log"
	    
        echo "input: $f1, output: $f2"			
	    /Datos/tools/usearch11.0.667_i86linux32 -threads 16 -fastx_uniques ./$f1 -fastaout ../UNIQUESbyLIB/$f2 -fastqout ../UNIQUESbyLIB/$f2_fq -sizeout &> ../UNIQUESbyLIB/$f2log
	     		
    done
    echo "done fastx_uniques of fastas ..."

    cat ../UNIQUESbyLIB/*.fas >> ../UNIQUESbyLIB/all_reads_UNIQUESbyLIB.fasta
    /Datos/tools/usearch11.0.667_i86linux32 -fastx_uniques ../UNIQUESbyLIB/all_reads_UNIQUESbyLIB.fasta -fastaout ../UNIQUESbyLIB/all_reads_UNIQUESbyLIB_derep.fasta -sizeout

    echo "executing fastx_uniques of fastqs..."
    for f1 in *sorted416_420.fq
    do
	    
	    #echo $f1
	    part1="$(echo $f1 | sed 's/\(.*\)fq*/\1/')" 

        f2_fq="$part1"
	    f2_fq+="uniques.fq"
        f2log="$part1"
	    f2log+="uniques.log"
	    
        echo "input: $f1, output: $f2_fq"		

	    /Datos/tools/usearch11.0.667_i86linux32 -threads 16 -fastx_uniques ./$f1 -fastqout ../UNIQUESbyLIB/$f2_fq -sizeout

    done			
    echo "done fastx_uniques of fastqs..."
fi


##############################################################################################################################
#                                                                                                                       
#                                           UNOISE 3
#
##############################################################################################################################


# unoise3 by lib
# error correction
# denoised zero-radius OTUs (zOTUs)
# minsize 2 to include all sequences minus singletons
# tabbedout to identifiy chimeras and eliminate them

# test different alpha values for unoise
if [ ! -d "../../03_denoising/UNOISE_a1" ]
then
    cd ../UNIQUESbyLIB
    for i in {1..10}
    do
        if [ ${i} != 2 ]
        then 
            mkdir ../../03_denoising/UNOISE_a${i}
            alpha="${i}.0"
            echo "executing unoise3 with alpha ${alpha}..."

            for f1 in *uniques.fas
            do
                
                #echo $f1
                part1="$(echo $f1 | sed 's/\(.*\)fas*/\1/')" 			
                f2="$part1"
                f2+="zotus.fas"
                f2log="$part1"
                f2log+="zotus.log"
                f2tabout="$part1"
                f2tabout+="zotus.tabbedout.txt"
                                    
                echo "input: $f1, output: $f2"

                /Datos/tools/usearch11.0.667_i86linux32 -unoise3 ./$f1 -zotus ../../03_denoising/UNOISE_a$i/$f2 -minsize 2 -threads 16 -unoise_alpha $alpha -tabbedout ../../03_denoising/UNOISE_a$i/$f2tabout &>../../03_denoising/UNOISE_a$i/$f2log
                
            done
            echo "done unoise3 ..."
        fi
    done
fi


# RELABEL and join UNOISE
if [ ! -d "../../03_denoising/UNOISE_relabelled_joined_a1" ]
then
    echo "Executing UNOISE with different alpha values.."
    for i in {1..10}
    do
        if [ ${i} != 2 ]
        then 
            cd ../../03_denoising/UNOISE_a$i
            mkdir ../UNOISE_relabelled_joined_a$i
            for file in *zotus.fas
            do	
                part1="$(echo $file | sed 's/\(.*\)_R1.*.*/\1/')" 
                label="$part1"
                
                echo "input: $file, output: $label"
                sed "-es/^>\(.*\)/>\1_$label/" < ${file} > ..//UNOISE_relabelled_joined_a$i/${file%.*}.relabel.fas
            
            done

            # Join all zotus individuals
            cd ../UNOISE_relabelled_joined_a$i
            cat *relabel.fas >> all.ZOTUSbySAMPLE.fas
            rm *relabel.fas

            # Only keep uniques (in the fasta that inlcudes all inviduals)
            /Datos/tools/usearch11.0.667_i86linux32 -fastx_uniques all.ZOTUSbySAMPLE.fas -fastaout all.ZOTUS_UNIQUES.fas -threads 16 &> all.ZOTUS_UNIQUES.log

            vsearch --search_exact ../../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db all.ZOTUS_UNIQUES.fas --sizeout --dbmatched all.ZOTUS_UNIQUES_with_size.fasta --threads 10
            vsearch --uchime3_denovo all.ZOTUS_UNIQUES_with_size.fasta -uchimeout uchime_out.txt -chimeras uchime_chimeras.fa -nonchimeras unoise3_uchime.fasta -xsize

            kraken2 --db /Datos/tools/kraken2_classification/kraken2_db/kraken2_nt/ --threads 10 \
            --output unoise_kraken2_nt.txt unoise3_uchime.fasta --memory-mapping --use-names --report unoise_kraken2_nt_report.txt

            python3 /Datos/tools/KrakenTools/extract_kraken_reads.py -k unoise_kraken2_nt.txt \
            -s unoise3_uchime.fasta -o unoise3_uchime_onlyArthropods.fasta -t 6656 --include-children -r unoise_kraken2_nt_report.txt

            vsearch --search_exact ../../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db unoise3_uchime_onlyArthropods.fasta --otutabout ASV_table_unoise3.csv --threads 10
            vsearch --usearch_global unoise3_uchime_onlyArthropods.fasta --db unoise3_uchime_onlyArthropods.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads 10

        fi
    done
fi


##############################################################################################################################
#
#                                               dnoise (2021)
#
##############################################################################################################################
cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/UNIQUESbyLIB


if [ ! -d "../preprocess_fasta_for_dnoise" ]
then
    echo "generating files for dnoise.."
    mkdir ../preprocess_fasta_for_dnoise
    cd ../preprocess_fasta_for_dnoise
    cp  ../UNIQUESbyLIB/*.fas .

    # dnoise does not accept multiple semicolons 
    # there needs to be a semicolon at the end at to separate size from the rest of the header    
    find *.fas -exec bash -c "sed -i 's/;b/_b/g' {}" \;

fi

# first nt codon position is -x 3
# -y entropy
# -g to only return entropy levels from dataset
# -j joining criteria
# --alpha default 5
# --min_abund (default 1)
# --first_nt_codon_position (-x)
# --count_name ('size' by default)
cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/preprocess_fasta_for_dnoise

# testing different alpha values for dnoise
cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/preprocess_fasta_for_dnoise
if [ ! -d "../../03_denoising/dnoise_alpha1" ]
then
    for i in {1..10}
    do
        mkdir ../../03_denoising/dnoise_alpha${i}
        echo "Executing dnoise with alpha value $i.."

        for f1 in *uniques.fas
        do
            echo $f1
            dnoise -c 18 --fasta_input $f1 --fasta_output ../../03_denoising/dnoise_alpha${i}/$f1"_dnoise_default_alpha${1}_y.fasta" -y --first_nt_codon_position 3 --min_abund 2 --alpha ${i}
        done
    done
fi

# RELABEL and join dnoise
if [ ! -d "../../03_denoising/dnoise_relabelled_joined_alpha1" ]
then
    conda activate kraken2_env
    echo "Executing dnoise relabelling, chimera removal, non-arthopod removal.."
    for i in {1..10}
    do
        cd ../../03_denoising/dnoise_alpha${i}
        mkdir ../../03_denoising/dnoise_relabelled_joined_alpha${i}

        # Join all zotus individuals
        cat *.fasta >> ../dnoise_relabelled_joined_alpha${i}/all.MOTUSbySAMPLE.fasta
        cd ../dnoise_relabelled_joined_alpha${i}

        # Only keep uniques (in the fasta that inlcudes all inviduals)
        /Datos/tools/usearch11.0.667_i86linux32 -fastx_uniques all.MOTUSbySAMPLE.fasta -fastaout all.MOTUS_UNIQUES.fasta -threads 10 -relabel ASV &> all.MOTUS_UNIQUES.log
        vsearch --search_exact ../../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db all.MOTUS_UNIQUES.fasta --sizeout --dbmatched all.MOTUS_UNIQUES_with_size.fasta --threads 10
        vsearch --uchime3_denovo all.MOTUS_UNIQUES_with_size.fasta -uchimeout uchime_out.txt -chimeras uchime_chimeras.fa -nonchimeras dnoise_uchime.fasta -xsize

        kraken2 --db /Datos/tools/kraken2_classification/kraken2_db/kraken2_nt/ --threads 10 \
        --output dnoise_kraken2_nt.txt dnoise_uchime.fasta --memory-mapping --use-names --report dnoise_kraken2_nt_report.txt

        python3 /Datos/tools/KrakenTools/extract_kraken_reads.py -k dnoise_kraken2_nt.txt \
        -s dnoise_uchime.fasta -o dnoise_uchime_onlyArthropods.fasta -t 6656 --include-children -r dnoise_kraken2_nt_report.txt

        vsearch --search_exact ../../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db dnoise_uchime_onlyArthropods.fasta --otutabout ASV_table_dnoise.csv --threads 10
        vsearch --usearch_global dnoise_uchime_onlyArthropods.fasta --db dnoise_uchime_onlyArthropods.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads 10
    done
fi




##############################################################################################################################
#
#                                               DADA2 (2016)
#
##############################################################################################################################
# --p-trunc-len 0 = no truncation will be performed

if [ ! -d "../dada2_omega1e-90" ]
then
    echo "Running dada2.."
    Rscript --vanilla /Datos/tools/dada2_pip.R /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/dada2_omega1e-90/ 1e-90 19
    vsearch --usearch_global dada2_default.fasta --db dada2_default.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads 10  
    vsearch --search_exact ../../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db dada2_default.fasta --otutabout ASV_table_dada2.csv --sizeout --dbmatched dada2_default_with_size.fasta --threads 10
    vsearch --uchime3_denovo all.MOTUS_UNIQUES_with_size.fasta -uchimeout uchime_out.txt -chimeras uchime_chimeras.fa -nonchimeras all.MOTUS_UNIQUES_with_size_uchime.fasta -xsize

fi



if [ ! -d "../dada2_omega1e-90" ]
then

    cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/03_denoising/


    echo "Executing DADA2 chimera removal .."

    for D in $(find dada2_omega*/dada2_*.fasta | grep -v 'chimera');
    do
        parentPath="$(dirname "$D")"
        parentDir="$(basename $parentPath)"
        echo "Running kraken2 on dada2.."
        vsearch --search_exact ../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db $D --sizeout --dbmatched $parentDir/dada2_with_size.fasta --threads 10
        vsearch --uchime3_denovo $parentDir/dada2_with_size.fasta -uchimeout $parentDir/uchime_out.txt -chimeras $parentDir/uchime_chimeras.fa -nonchimeras $parentDir/dada2_uchime.fasta -xsize
        
        kraken2 --db /Datos/tools/kraken2_classification/kraken2_db/kraken2_nt/ --threads 18 \
        --output $parentDir/dada2_lulu_kraken2_nt.txt $D --memory-mapping --use-names --report $parentDir/dada2_lulu_kraken2_nt_report.txt

        python3 /Datos/tools/KrakenTools/extract_kraken_reads.py -k $parentDir/dada2_lulu_kraken2_nt.txt \
        -s $D -o $parentDir/dada2_onlyArthropods.fasta -t 6656 --include-children -r $parentDir/dada2_lulu_kraken2_nt_report.txt

        vsearch --search_exact ../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db $parentDir/dada2_onlyArthropods.fasta --otutabout $parentDir/ASV_table_dada2.csv --threads 10
        vsearch --usearch_global $parentDir/dada2_onlyArthropods.fasta --db $parentDir/dada2_onlyArthropods.fasta --self --id .84 --iddef 1 --userout $parentDir/match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads 10
    done

fi

##############################################################################################################################
#
#                                               SWARM 
#
##############################################################################################################################
# execute swarm 
# swarm -d default = 1
cd /Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/UNIQUESbyLIB
if [ ! -d "../../03_denoising/swarm_d1" ]
then
    echo "Running swarm.."

    for i in {1..10}
    do
        mkdir ../../03_denoising/swarm_d${i}
        echo "Executing swarm with d = ${i}.."

        for f1 in *uniques.fas
        do
            echo $f1
            sed '/^>/ ! s/[Nn]//g' $f1 | swarm -t 18 -d ${i} -o ../../03_denoising/swarm_d${i}/$f1"_swarm_output.txt" -u ../../03_denoising/swarm_d${i}/$f1"_swarm_uclust_format.txt" -w ../../03_denoising/swarm_d${i}/$f1"_swarm_representatives.fasta" -z
        done
    done

fi

# relabel swarm
if [ ! -d "../../03_denoising/swarm_relabelled_d1" ]
then
    echo "Relabelling SWARM .."
    for i in {1..10}
    do
        mkdir ../../03_denoising/swarm_relabelled_d${i}
        cd ../../03_denoising/swarm_d${i}

        # Join all cluster representatives
        cat *.fasta >> ../swarm_relabelled_d${i}/all.clusterRepresentativesBySAMPLE.fasta
        cd ../swarm_relabelled_d${i}

        # Only keep uniques (in the fasta that inlcudes all inviduals)
        /Datos/tools/usearch11.0.667_i86linux32 -fastx_uniques all.clusterRepresentativesBySAMPLE.fasta -fastaout all.cluster_UNIQUES.fasta -threads 10 &> all.cluster_UNIQUES.log
        vsearch --sortbysize all.cluster_UNIQUES.fasta --minsize 2 --output all.cluster_UNIQUES_minAbund2.fasta
        sed -E 's/;barcodelabel=[^;]*;size=[0-9]+;//g' all.cluster_UNIQUES_minAbund2.fasta > all.cluster_UNIQUES_minAbund2_renamed.fasta
        vsearch --search_exact ../../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db all.cluster_UNIQUES_minAbund2_renamed.fasta --sizeout --dbmatched swarm_minAbund2_with_size.fasta --threads 10
        vsearch --uchime3_denovo swarm_minAbund2_with_size.fasta -uchimeout uchime_out.txt -chimeras uchime_chimeras.fa -nonchimeras swarm_minAbund2_uchime.fasta -xsize
        vsearch --usearch_global swarm_minAbund2_uchime.fasta --db swarm_minAbund2_uchime.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads 10
        # sed -E 's/;size=[0-9]+;//g' match_list.txt > match_list_edited.txt
        vsearch --search_exact ../../02_preprocessed/LENGTH416_420/all_reads_LENGTH416_420.fasta  --db swarm_minAbund2_uchime.fasta --otutabout ASV_table_swarm.csv --threads 10

    done
fi
