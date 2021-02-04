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
known_seq_path = os.path.join('resources', 'known_seqs')

# Locally, the known-sequence set fits in memory. On the cloud, it's potentially terabytes, and is built and run totally differently.
running_on_cloud = '--cloud' in sys.argv
if running_on_cloud:
    if not running_as_script:
        raise NotImplementedError('the --cloud option is currently only designed to work when running as a script.')
    import multiprocessing
    cpus = multiprocessing.cpu_count()
else:
    # create the known_seq_set as an importable module-level variable (not cloud-scale)
    known_seq_set = set()
    print('building known_seq fragments to exclude...')
    for dir_entry in os.listdir(known_seq_path):
        if tuple(dir_entry.split('.')[-2:]) != ('fragset', 'json'):
            continue
        with open(os.path.join(known_seq_path, dir_entry)) as f:
            known_seq_set.update(json.loads(f.read()))
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
                if not running_on_cloud and variant in known_seq_set:
                    print('not including', variant)
                    continue
                hazard_set.append(variant)
    print('done')

    if running_on_cloud:

        def split_into_windows(str_in, window=19):
            for i in range(len(str_in)):
                if i + window <= len(str_in):
                    yield str_in[i:i+window]

        approx_batch_size = 8e6 # a few dozen megabytes at a time
        hazard_set = frozenset(hazard_set)
        with open('.' + output_path + '_temp', 'w+') as f:
            f.write(json.dumps(hazard_set, indent=2)) #TODO: do we want to be generating this temp file every time?

        def known_seqs_in_hazard_set(known_seq_fname):
            from Bio import SeqIO
            window_batch = set()
            intersect_set = set()
            for record in SeqIO.parse(known_seq_fname, "fasta"): # potentially many GBs
                if len(window_batch) < approx_batch_size:
                    window_batch.update(split_into_windows(str(record.seq)))
                    continue
                intersect_set.update(set.intersection(hazard_set, window_batch))
                window_batch = set()
            intersect_set.update(set.intersection(hazard_set, window_batch))
            return intersect_set

        with multiprocessing.Pool(processes=cpus) as pool:
            known_hazard_seqs = set.union(pool.imap_unordered(known_seqs_in_hazard_set, os.listdir(known_seq_path)))
        hazard_set = hazard_set - known_hazard_seqs

    print('writing output...')
    with open(output_path, 'w+') as f:
        f.write(json.dumps(hazard_set, indent=2))
    print('done')

