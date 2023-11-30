#!/bin/bash/

#Bash script, maps collapsed smallRNA-seq fasta files to MirGeneDB 
#extended mature/star bowtie reference. Uses bowtie 1.2, allowing up
#to 10 valid alignments, from best to worst, using -v 3 mismatches 
#policy, of allowing up to 3 mismatches in the alignment. Outputs 
#sam files.

for k in /home/jcdenton/projects/smallRNAseq/*

do
    echo ${k##*/}
    rm -r /home/jcdenton/projects/smallRNAseq_sam/${k##*/}
    mkdir /home/jcdenton/projects/smallRNAseq_sam/${k##*/}
    > /home/jcdenton/projects/smallRNAseq_sam/${k##*/}/mirgenedb_${k##*/}.log
    for i in $k/*.fa*
    do
        echo ${i##*/}
        bowtie -f -k 10 --threads 4 --best --norc -v 3 /home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/${k##*/}.fas \
        $k/${i##*/} -S > /home/jcdenton/projects/smallRNAseq_sam/${k##*/}/${i##*/}.sam 2>> /home/jcdenton/projects/smallRNAseq_sam/${k##*/}/mirgenedb_${k##*/}.log
    done
done
