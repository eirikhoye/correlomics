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
