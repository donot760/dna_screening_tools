import json
import os
import sys
from Bio.Seq import Seq
from ecoli_seq_variants import single_variants
from metro_hastings_variants import MetroHastingsVariants
from funtrp_blosum_predict import parse_funtrp_line

running_as_script = __name__ == '__main__'
if running_as_script and len(sys.argv) < 2:
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

if running_as_script:
    hazard_set = []
    print('generating and filtering variants...')
    with open('resources/aa_fragment_picks.txt') as f:
        for funtrp_line in f:
            funtrp_line = funtrp_line.strip()
            if not funtrp_line:
                continue
            aa_frag, ntr_triples = parse_funtrp_line(funtrp_line)
            single_variant_set = set((str(Seq(v).translate()) for v in single_variants(aa_frag)))
            #print('num single variants:', len(single_variant_set))
            variant_computer = MetroHastingsVariants(aa_frag, ntr_triples, 10000)
            all_variants = single_variant_set | variant_computer.result_set()
            #print('num Metro-Hastings variants:', len(variant_computer.result_set()))
            for variant in all_variants: # includes original (zero-variant) too
                variant = str(variant)
                if variant in genbank_set:
                    print('not including', variant)
                    continue
                hazard_set.append(variant)
    print('done')

    print('writing output...')
    with open(output_path, 'w+') as f:
        f.write(json.dumps(hazard_set, indent=2))
    print('done')

