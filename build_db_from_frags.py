import json
import os
import sys
import time
from Bio.Seq import Seq
from ecoli_seq_variants import single_variants
from metro_hastings_variants import MetroHastingsVariants
from funtrp_blosum_predict import parse_funtrp_line

running_as_script = __name__ == '__main__'
known_seq_path = os.path.join('resources', 'known_seqs')
if running_as_script and len(sys.argv) < 2:
    print('must supply output file path as last argument')
    print('options:\n\t--cloud: run with thread pool and possibly large db'
          '\n\t--known-seq-path: path to find database of known sequences to be'
          ' pre-screened out of the hazard variants database. default: '
          + known_seq_path)
    exit()

output_path = sys.argv[-1]
try:
    known_seq_path = sys.argv[sys.argv.index('--known-seq-path') + 1]
except ValueError:
    pass

start_time = time.time()

# Locally, the known-sequence set fits in memory. On the cloud, it's potentially terabytes, and is built and run totally differently.
running_on_cloud = '--cloud' in sys.argv
if running_on_cloud:
    if not running_as_script:
        raise NotImplementedError('the --cloud option is currently only designed to work when running as a script.')
    import multiprocessing
    cpus = multiprocessing.cpu_count()
    print('using', cpus, 'cores.')
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

original_frags = set()

if running_as_script:
    hazard_set = []
    print('generating and filtering variants...')
    if os.path.isfile('resources/hazard_scinames.txt'):
        with open('resources/hazard_scinames.txt') as f:
            hazard_scinames = [line.strip().lower() for line in f]
    else:
        hazard_scinames = []
    with open('resources/aa_fragment_picks.txt') as f:
        for funtrp_line in f:
            funtrp_line = funtrp_line.strip()
            if not funtrp_line:
                continue
            aa_frag, ntr_triples = parse_funtrp_line(funtrp_line)
            original_frags.add(aa_frag)
            single_variant_set = set((str(Seq(v).translate()) for v in single_variants(aa_frag)))
            #print('num single variants:', len(single_variant_set))
            variant_computer = MetroHastingsVariants(aa_frag, ntr_triples, 1000000)
            all_variants = single_variant_set | variant_computer.result_set()
            #print('num Metro-Hastings variants:', len(variant_computer.result_set()))
            for variant in all_variants: # includes original (zero-variant) too
                variant = str(variant)
                if not running_on_cloud and variant in known_seq_set:
                    print('not including', variant)
                    continue
                hazard_set.append(variant)
    print(time.time() - start_time)
    tot_m, used_m, free_m = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
    print(used_m/tot_m*100, '% memory used')
    print('done')

    if running_on_cloud:

        def split_into_windows(str_in):
            if not set(str_in.lower()) - set('atcg \t\n'): # if it's DNA:
                if False: # enable in the future when actually screening DNA
                    for i in range(len(str_in)):
                        if i + 42 <= len(str_in):
                            yield str_in[i:i+42]
                seq = Seq(str_in)
                for sense in seq, seq.reverse_complement():
                    for frame_offset in range(3):
                        frame = str(sense[frame_offset:].translate())
                        for i in range(len(frame)):
                            if i + 19 <= len(frame):
                                yield frame[i:i+19]
            else: # it's a protein/translation
                for i in range(len(str_in)):
                    if i + 19 <= len(str_in):
                        yield str_in[i:i+19]

        approx_batch_size = 8 # a few dozen megabytes at a time
        hazard_set = frozenset(hazard_set)
        with open(os.path.join(os.path.dirname(output_path), '.' + os.path.basename(output_path) + '_temp'), 'w+') as f:
            f.write(json.dumps(list(hazard_set), indent=2)) #TODO: do we want to be generating this temp file every time?

        def known_seqs_in_hazard_set(known_seq_fname):
            #from Bio import SeqIO
            #window_batch = set()
            intersect_set = set()
            seq_list = []
            #for record in SeqIO.parse(known_seq_fname, "fasta"): # potentially many GBs
            skipping = False
            with open(known_seq_fname) as f: # potentially many GBs
                for line in f:
                    if not line:
                        continue
                    if line[0] == '>':
                        seq = ''.join(seq_list)
                        for window in split_into_windows(seq):
                            if window in hazard_set:
                                intersect_set.add(window)
                        line = line.lower()
                        skipping = any(sciname in line for sciname in hazard_scinames) # go to skip mode to skip over the following sequence
                        seq_list = []
                        continue
                    if skipping:
                        continue
                    seq_list.extend(line.strip())
            # don't forget leftover/partially complete items
            seq = ''.join(seq_list)
            for window in split_into_windows(seq):
                if window in hazard_set:
                    print(window, '<- that window was present')
                    intersect_set.add(window)
            if intersect_set:
                print(known_seq_fname + ':', intersect_set)
            return intersect_set

        indiv_paths = [os.path.join(known_seq_path, fname) for fname in os.listdir(known_seq_path) if fname[0] != '.']
        indiv_paths = [p for p in indiv_paths if os.path.isfile(p)]
        print(indiv_paths)

        print('initial size of hazard set:', len(hazard_set))
        with multiprocessing.Pool(processes=cpus) as pool:
            known_hazard_seqs = set.union(*pool.imap_unordered(known_seqs_in_hazard_set, indiv_paths))
        print('removing', len(known_hazard_seqs), 'known hazard sequence fragment(s)')
        if len(known_hazard_seqs) < 100:
            print(known_hazard_seqs)
        hazard_set = (hazard_set - known_hazard_seqs) | original_frags
        print('final size of hazard set:', len(hazard_set))

    print('writing output...')
    with open(output_path, 'w+') as f:
        f.write(json.dumps(list(hazard_set), indent=2))
    print('done')

