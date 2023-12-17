Correlomics:uncovering microRNA editing and structural features based on reference and sequencing data

Pipeline for smallRNA sequencing data. 

Processing and quality control of sequencing data using MirTrace. Returns fasta files with collapsed reads.

Generates microRNA reference library from MirGeneDB mature and star microRNA sequences, which are extended with 5 nucleotides at the 5p and 3p end from the genomic loci. 

Generate bowtie reference for each species, bowtie v=1.0.0. 

Map collapsed reads, allowing up to 10 alignments with up to 3 mismatches. Outputs sam files.

Use two separate filtering methods for sam files, to filter for multiple read alignments:
- One selecting for mismatches INSIDE miRNA position 2-18
    - Used for identifying microRNA editing events, including A to I editing.
- One selecting for mismatches OUTSIDE miRNA position 2-18
    - Used for analyzing the preferred read alignment position (0 is canonical start/end).


![Image](figures/flow_chart_editing_events.drawio.png) 





TODO
- Finish snakemake integration with mfold script
- Merge all reduntant scripts
- Update pipeline diagram
- Test on latest version of MirGeneDB with smallRNA-seq at SAGA