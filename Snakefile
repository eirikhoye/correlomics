rule extend_mature_star_pre_ref:
    """
    Extend all MirGeneDB mature and star sequences with 5nt 5p and 3p, 
    using precursors as reference
    """
    input:
        