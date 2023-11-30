# BASH script that builds a bowtie reference from extended MirGeneDB mature  star sequences

for k in /home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_genome/*

do 
    echo $k
    bowtie-build $k /home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/${k##*/}
done

