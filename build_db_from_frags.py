import json
import os
import sys
from Bio.Seq import Seq
from ecoli_seq_variants import single_variants

if len(sys.argv) < 2:
    print('must supply output file path')
    exit()

output_path = sys.argv[-1]

genbank_path = os.path.join('resources', 'genbank')
genbank_set = set()
print('building genbank fragments to exclude...')
for dir_entry in os.listdir(genbank_path):
    if tuple(dir_entry.split('.')[-2:]) != ('fragset', 'json'):
        continue
    with open(os.path.join(genbank_path, dir_entry)) as f:
        genbank_set.update(json.loads(f.read()))
print('done')

hazard_set = []
print('generating and filtering variants...')
with open('resources/aa_fragment_picks.txt') as f:
    for aa_frag in f:
        aa_frag = aa_frag.strip()
        if not aa_frag:
            continue
        for variant in single_variants(aa_frag): # includes original (zero-variant) too
            variant_aa_frag = str(Seq(variant).translate())
            if variant_aa_frag in genbank_set:
                print('not including', variant_aa_frag)
                continue
            hazard_set.append(variant_aa_frag)
print('done')

print('writing output...')
with open(output_path, 'w+') as f:
    f.write(json.dumps(hazard_set))
print('done')

