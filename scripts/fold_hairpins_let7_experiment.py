import pickle
import os
import sys
from Bio import SeqIO

current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)

# Now you can import your module
from correlomics.foldings import Fold


"""
Define input files
"""
stemfile   = open('/home/jcdenton/projects/mirgenedb_database/merged_all-pri-30-30.fas', 'r')
maturefile = open('/home/jcdenton/projects/mirgenedb_database/mature.fas', 'r') 
starfile   = open('/home/jcdenton/projects/mirgenedb_database/star.fas')
loopfile   = open('/home/jcdenton/projects/mirgenedb_database/loopfile.fas')
stem   = SeqIO.to_dict(SeqIO.parse(stemfile, "fasta"))
mature = SeqIO.to_dict(SeqIO.parse(maturefile, "fasta"))
star   = SeqIO.to_dict(SeqIO.parse(starfile, "fasta"))
loop   = SeqIO.to_dict(SeqIO.parse(loopfile, "fasta"))


outdir = {
    'no_constraints': '/home/jcdenton/projects/correlomics/test_let7_no_constraints/',
    'restrict_lower_stem': '/home/jcdenton/projects/correlomics/test_let7_restrict_lower_stem/',
    'restrict_lower_stem_loop': '/home/jcdenton/projects/correlomics/test_let7_restrict_lower_stem_loop/',
    'restrict_loop': '/home/jcdenton/projects/correlomics/test_let7_restrict_loop/'
}


outdir = {
    'no_constraints': '/home/jcdenton/projects/correlomics/test_no_constraints/',
    'restrict_lower_stem': '/home/jcdenton/projects/correlomics/test_restrict_lower_stem/',
    'restrict_lower_stem_loop': '/home/jcdenton/projects/correlomics/test_restrict_lower_stem_loop/',
    'restrict_loop': '/home/jcdenton/projects/correlomics/test_restrict_loop/'
}


# Process hairpins
for id in stem.keys():
#   if "Hsa-Let-7" not in id:
#       continue

    if "Mir-451" in id:
        continue

    hairpin = Fold(id[:-4])
    hairpin.fold_hairpin()
    with open(outdir['no_constraints'] + hairpin.mirgenedb_id + '.pkl', 'wb') as f:
        pickle.dump([hairpin.mirgenedb_id+'\n', hairpin.seq, hairpin.dot, hairpin.dG], f)

    hairpin = Fold(id[:-4])
    hairpin.fold_hairpin(constrain_stem=True)
    with open(outdir['restrict_lower_stem'] + hairpin.mirgenedb_id + '.pkl', 'wb') as f:
        pickle.dump([hairpin.mirgenedb_id+'\n', hairpin.seq, hairpin.dot, hairpin.dG], f)

    hairpin = Fold(id[:-4])
    hairpin.fold_hairpin(constrain_stem=True, constrain_loop=True)
    with open(outdir['restrict_lower_stem_loop'] + hairpin.mirgenedb_id + '.pkl', 'wb') as f:
        pickle.dump([hairpin.mirgenedb_id+'\n', hairpin.seq, hairpin.dot, hairpin.dG], f)

    hairpin = Fold(id[:-4])
    hairpin.fold_hairpin(constrain_stem=False, constrain_loop=True)
    with open(outdir['restrict_loop'] + hairpin.mirgenedb_id + '.pkl', 'wb') as f:
        pickle.dump([hairpin.mirgenedb_id+'\n', hairpin.seq, hairpin.dot, hairpin.dG], f)
