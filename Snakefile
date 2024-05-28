configfile: "config.yaml"

rule extend_mature_star_pre_ref:
    """
    Extend all MirGeneDB mature and star sequences with 5nt 5p and 3p, 
    using precursors as reference
    """
    input:
        config['database'] + config['prifile'],
        config['database'] + config['ref'] + "{species_id}.fas"
    output:
        config['database'] + config['extend_ref'] + "{species_id}.fas"
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/reference_extended_mature_star_maker.py {input} {output}"

rule bowtie_reference_maker:
    """
    Build bowtie references from extended MirGeneDB mature / star sequences
    """
    input:
        config['database'] + config['extend_ref'] + "{species_id}.fas"
    output:
        config['database'] + config['bowtie_ref'] + "{species_id}.1.ebwt",
        config['database'] + config['bowtie_ref'] + "{species_id}.2.ebwt",
        config['database'] + config['bowtie_ref'] + "{species_id}.3.ebwt",
        config['database'] + config['bowtie_ref'] + "{species_id}.4.ebwt",
        config['database'] + config['bowtie_ref'] + "{species_id}.rev.1.ebwt",
        config['database'] + config['bowtie_ref'] + "{species_id}.rev.2.ebwt"
        
    conda:
        "envs/bowtie.yaml"
    shell:
        "bowtie-build {input} {config['database']}/mirgenedb_merged_extended_reference/{wildcards.species_id}"

rule bowtie_mirgenedb_collapsed_species_merged_extended_v2:
    """
    Maps collapsed smallRNA-seq fasta files to MirGeneDB extended mature/star bowtie reference. Uses
    bowtie 1.2, allowing up to 10 valid alignments, from best to worst, using -v 3 mismatches policy, of allowing up to 3 mismatches in the alignment. 
    Outputs sam files.
    """
    input:
        reference = config['database'] + config['bowtie_ref'] + "{species_id}.1.ebwt",
        fastafile = config['seq_data'] + "{species_id}/{tissue_id}.fas"
    params:
        config['database'] + config['bowtie_ref'] + "{species_id}"
    output:
        samfile = config['database'] + config['sam'] + "{species_id}/{tissue_id}.sam"
    log:
        config['database'] + config['sam'] + "log/{species_id}/{tissue_id}.log"
    conda:
        "envs/bowtie.yaml"
    shell:
        "bowtie -f -k 10 --best --norc -v3 {params} {input.fastafile} -S > {output.samfile} 2> {log}"

rule filter_samfile_no_internal_mismatch:
    input:
        config['database'] + config['sam'] + "{species_id}/{tissue_id}.sam"
    output:
        config['database'] + config["filtered_no_pos_2_18_mis"] + "{species_id}/{tissue_id}_filtered.sam"
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/filter.samfile.species.merged.extended.no.position.penalty_no_mismatch_pos2_18.py {input} {output}"

rule filter_samfile_only_internal_mismatch:
    input:
        config['database'] + config['sam'] + "{species_id}/{tissue_id}.sam"
    output:
        config['database'] + config["filtered_no_outside_pos_2_18_mis"] + "{species_id}/{tissue_id}_filtered.sam"
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/filter.samfile.species.merged.extended.no.position.penalty_no_mismatch_outside_pos2_18.py {input} {output}"

rule count_editing_events:
    input:
        config['species_dict'],
        config['database'] + config["extend_ref"] + "{species_id}.fas",   
        config['database'] + config["filtered_no_outside_pos_2_18_mis"] + "{species_id}/{tissue_id}_filtered.sam"
    output:
        config['database'] + config['editing_events'] + '{species_id}/{tissue_id}_editing_events.csv'
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/samfile.all.nucleotides.edit.counter.mirna.py {input} {output} {wildcards.species_id} {wildcards.tissue_id}"

rule aggregate_editing_events:
    input:
        config['database'] + config['editing_events'] + '{species_id}/'
    output:
        config['database'] + config['editing_events_aggregated'] + '{species_id}_editing_aggregated.csv'
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/concat_editing_events.py {input} {output}"

rule collapse_mismatches_species:
    input:
        config['database'] + config['editing_events_aggregated'] + '{species_id}_editing_aggregated.csv'
    output:
        config['database'] + config['counted_mismatches_collapsed_species'] + '{species_id}_editing_collapsed_species.csv'
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/collapse_mismatches.py {input} {output} species"


rule collapse_mismatches_tissue:
    input:
        config['database'] + config['editing_events_aggregated'] + '{species_id}_editing_aggregated.csv'
    output:
        config['database'] + config['counted_mismatches_collapsed_tissue'] + '{species_id}_editing_collapsed_tissue.csv'
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/collapse_mismatches.py {input} {output} tissue"
