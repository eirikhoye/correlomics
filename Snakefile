#configfile: "config.yaml"

rule extend_mature_star_pre_ref:
    """
    Extend all MirGeneDB mature and star sequences with 5nt 5p and 3p, 
    using precursors as reference
    """
    input:
        "/home/jcdenton/projects/mirgenedb_database/merged_all-pri-30-30.fas",
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_mature_star_genome/{species_id}.fas"
    output:
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_genome/{species_id}.fas"
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/reference_extended_mature_star_maker.py {input} {output}"

rule bowtie_reference_maker:
    """
    Build bowtie references from extended MirGeneDB mature / star sequences
    """
    input:
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_genome/{species_id}.fas"
    output:
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/{species_id}.1.ebwt",
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/{species_id}.2.ebwt",
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/{species_id}.3.ebwt",
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/{species_id}.4.ebwt",
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/{species_id}.rev.1.ebwt",
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/{species_id}.rev.2.ebwt"
        
    conda:
        "envs/bowtie.yaml"
    shell:
        "bowtie-build {input} /home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/{wildcards.species_id}"

rule bowtie_mirgenedb_collapsed_species_merged_extended_v2:
    """
    Maps collapsed smallRNA-seq fasta files to MirGeneDB extended mature/star bowtie reference. Uses
    bowtie 1.2, allowing up to 10 valid alignments, from best to worst, using -v 3 mismatches policy, of allowing up to 3 mismatches in the alignment. 
    Outputs sam files.
    """
    input:
        "/home/jcdenton/projects/smallRNAseq/{species_id}/{tissue_id}.fas"
    output:
        "/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_sam/{species_id}/{tissue_id}.sam"
    conda:
        "envs/bowtie.yaml"
    shell:
        "bowtie -f -k 10 --best --norc -v3 /home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_reference/{wildcards.species_id} {input} -S > {output} 2>> /home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_sam/{wildcards.species_id}/mirgenedb_{wildcards.species_id}.log"