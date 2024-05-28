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
stemfile   = open('/home/jcdenton/projects/mirgenedb_database/precursor.fas', 'r')
maturefile = open('/home/jcdenton/projects/mirgenedb_database/mature.fas', 'r') 
starfile   = open('/home/jcdenton/projects/mirgenedb_database/star.fas')
loopfile   = open('/home/jcdenton/projects/mirgenedb_database/loopfile.fas')
stem   = SeqIO.to_dict(SeqIO.parse(stemfile, "fasta"))
mature = SeqIO.to_dict(SeqIO.parse(maturefile, "fasta"))
star   = SeqIO.to_dict(SeqIO.parse(starfile, "fasta"))
loop   = SeqIO.to_dict(SeqIO.parse(loopfile, "fasta"))
outdir = '/home/jcdenton/projects/mirgenedb_database/hairpins_let7_pre/'


### debug
outdir = '/home/jcdenton/projects/correlomics/test_let7_pre/'

###

# Open output files
output_files = {
    'no_constraints': open('misfold_no_constraints.txt', 'w'),
    'constrain_lower_stem': open('misfold_constrain_lower_stem.txt', 'w'),
    'constrain_lower_stem_and_loop': open('misfold_constrain_lower_stem_and_loop.txt', 'w')
}


def process_hairpin(hairpin, output_file):
    if hairpin.check_misfold() is False:
        with open(outdir + hairpin.mirgenedb_id + '.pkl', 'wb') as f:
            pickle.dump([hairpin.seq, hairpin.dot, hairpin.dG], f)
    else:
        output_file.write(hairpin.id + '\n')
        output_file.write(hairpin.seq)
        output_file.write(hairpin.dot)


# Process hairpins
for id in stem.keys():
    if "Mir-451" in id:
        continue

    if "Hsa-Let-7" not in id:
        continue



    hairpin = Fold(id[:-4])  # drop miRNA without loop seq
    if hairpin.mirgenedb_id + '_loop' not in loop.keys():
        continue

    hairpin.fold_hairpin()  # no constraints
    if not hairpin.check_misfold():
        process_hairpin(hairpin, output_files['no_constraints'])
    else:
        hairpin.fold_hairpin(constrain_stem=True)  # Constrain lower stem to not fold
        if not hairpin.check_misfold():
            process_hairpin(hairpin, output_files['constrain_lower_stem'])
        else:
            hairpin.fold_hairpin(constrain_stem=True, constrain_loop=True)  # Constrain lower stem and loop to not fold
            if not hairpin.check_misfold():
                process_hairpin(hairpin, output_files['constrain_lower_stem_and_loop'])

# Close output files
for file in output_files.values():
    file.close()