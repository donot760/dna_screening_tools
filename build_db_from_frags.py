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
          + known_seq_path +
          '\n\t--role: when using delegation to address memory concerns,'
          '\n\t\t"helper": save intermediate pre-screen results in parent dir'
          ' of known-seq-path. ignores output path.'
          '\n\t\t"aggregate": assuming helpers are finished, assemble final'
          ' result from intermediates and save to output path. Note: typically,'
          ' the helpers will work in subdirectories of the known-seq-path of'
          ' the aggregator.'
          '\n\t\tany other value or missing: delegation disabled')
    exit()

# last terminal argument becomes output path
output_path = sys.argv[-1]

try:
    # example: in $ python build_db_from_frags.py --known-seq-path ~/blastdb resources/hazard_db/hazards.fragset.json
    # "~/blastdb" is the known sequence database path
    known_seq_path = sys.argv[sys.argv.index('--known-seq-path') + 1]
except ValueError:
    pass

# Locally, the known-sequence set fits in memory. On the cloud, it's potentially terabytes, and is built and run totally differently.
# use the --cloud switch to run with parallelism (multiple cores) and a large database
running_on_cloud = '--cloud' in sys.argv

role = None
try:
    # example: in $ python build_db_from_frags.py --role helper resources/hazard_db/hazards.fragset.json
    # "helper" is the role name. the helper create intermediate results that are used by the "aggregate" role later.
    role = sys.argv[sys.argv.index('--role') + 1]
except ValueError:
    pass
if role is not None and not running_on_cloud:
    print('warn: roles are only used with --cloud. did you mean to use that option?')

# num_vars is the number of variants we will generate for each fragment in aa_fragment_picks.txt.
num_vars = 100000 # default value
try:
    num_vars = sys.argv[sys.argv.index('--num-vars') + 1]
except ValueError:
    pass
num_vars = int(num_vars)
print('number of variants generated for each aa fragment:', num_vars)

# cache_path specifies a file where the results of the metropolis-hastings algorithm (variant set) will be stored.
# This is in case we have to cancel later, so we don't lose the results of this intensive computation, OR if
# the helper role is in use, which reads this cache.
cache_path = os.path.join(os.path.dirname(output_path), '.' + os.path.basename(
        output_path) + '_' + str(num_vars) + 'vars_temp')

start_time = time.time() # only for debug & metrics

if running_on_cloud:
    if not running_as_script:
        raise NotImplementedError('the --cloud option is currently only designed to work when running as a script.')
    import multiprocessing # use multiple cores & parallelism in computations
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

# original_frags is the small set of fragments that appear literally in aa_fragment_picks.txt
# we'll add them back in at the end, in case they were pre-screened out.
original_frags = set()

if running_as_script:
    print('generating and filtering variants...')
    # a list of the hazards we are screening (the organisms from which the entries in aa_fragment_picks.txt are drawn)
    # are listed in resources/hazard_scinames.txt. this helps us not pre-screen out fragments of the actual hazards we're screening.
    if os.path.isfile('resources/hazard_scinames.txt'):
        with open('resources/hazard_scinames.txt') as f:
            # make sure the comparison is case-insensitive by always comparing in lowercase (lower())
            hazard_scinames = [line.strip().lower() for line in f]
    else:
        print('warning: no resources/hazard_scinames.txt found.'
              ' hazards will not be skipped based on their scientific names in'
              ' sequence annotations.')
        hazard_scinames = []
    aa_frags, ntrs = [], []
    with open('resources/aa_fragment_picks.txt') as f:
        for funtrp_line in f:
            funtrp_line = funtrp_line.strip()
            if not funtrp_line:
                continue
            # ntr_triples is a list of 3-tuples of probabilities (floats). (neutral, toggle, rheostat). each tuple sums to 1.0.
            aa_frag, ntr_triples = parse_funtrp_line(funtrp_line)
            original_frags.add(aa_frag)
            # TODO: do we need both original_frags and aa_frags?
            aa_frags.append(aa_frag)
            ntrs.append(ntr_triples)

    # hazard_set_for_frag will be evaluated on different arguments in multiple processes in parallel
    # must have just 1 argument for parallel pool function evaluation with pool.imap_unordered(), below
    def hazard_set_for_frag(funtup): 
        aa_frag, ntr_triples = funtup
        hazard_vars = set()
        # at a minimum, we add the "obvious" variants resulting from changing just one amino acid at a time.
        # single_variants gives DNA strings {a,t,c,g} corresponding to all possible replacements at one amino acid position at a time
        # TODO: could optimize by using amino acids instead of going to DNA as an intermediate & translating
        single_variant_set = set((str(Seq(v).translate()) for v in single_variants(aa_frag))) # includes original fragment (zero-variant) too
        #print('num single variants:', len(single_variant_set))
        # we will also add many "less obvious" variants that we think are likely to work based on FUNTRP predictions and BLOSUM62.
        variant_computer = MetroHastingsVariants(aa_frag, ntr_triples, num_vars)
        # most heavy computation occurs in the following line.
        # combine the obvious and less obvious variants.
        all_variants = single_variant_set | variant_computer.result_set()
        #print('num Metro-Hastings variants:', len(variant_computer.result_set()))
        for variant in all_variants: # still includes original fragment (zero-variant)
            variant = str(variant)
            if not running_on_cloud and variant in known_seq_set:
                print('not including', variant)
                continue
            hazard_vars.add(variant)
        return hazard_vars

    if running_on_cloud:
        if os.path.isfile(cache_path):
            # if there is already a variant set computed in the past, load it and don't rerun metropolis-hastings
            print('using cached hazard variant set at', cache_path)
            with open(cache_path) as cache_f:
                hazard_set = set(json.loads(cache_f.read()))
                print('loaded', len(hazard_set), 'cached hazard variants')
        else:
            print('no cache; recalculating variants')
            with multiprocessing.Pool(processes=cpus) as pool:
                hazard_set = set.union(*pool.imap_unordered(hazard_set_for_frag,
                        zip(aa_frags, ntrs)))
            with open(cache_path, 'w+') as f:
                f.write(json.dumps(list(hazard_set), indent=2))
    else:
        hazard_set = []
        for funtup in zip(aa_frags, ntrs):
            hazard_set.extend(hazard_set_for_frag(funtup))
        hazard_set = set(hazard_set)

    print('done. time since start (s):', time.time() - start_time)
    # for debugging purposes, print out how much system memory has been consumed.
    tot_m, used_m, free_m = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
    print(int(used_m/tot_m*10000)/100, '% memory used')
    print('done')

    if running_on_cloud:

        def split_into_windows(str_in):
            # str_in could be DNA, RNA or amino acid sequences from e.g. genbank, nucleotide, NR
            # make case-insensitive by converting to lowercase. if it's RNA, convert it to equivalent DNA.
            # Note: this is only ok because U is not a valid amino acid abbreviation.
            str_in = str_in.lower().replace('u', 't')
            if not set(str_in[:1000]) - set('atcgn \t\n'): # if it's DNA:
                if False: # enable in the future when actually screening DNA
                    for i in range(len(str_in)):
                        # SecureDNA DNA fragments are always 42 base pairs long
                        if i + 42 <= len(str_in):
                            yield str_in[i:i+42]
                # convert the input string, which we now know is DNA, to a Seq object from the biopython module
                seq = Seq(str_in)
                # "sense" as in positive-sense and negative-sense DNA: whether you're reading it forward or backward. we will consider both.
                for sense in seq, seq.reverse_complement():
                    # "frame" in the sense of an open reading frame. there are 3 options, because each codon is 3 DNA bases long.
                    for frame_offset in range(3):
                        frame = str(sense[frame_offset:].translate())
                        for i in range(len(frame)):
                            # SecureDNA protein fragments are always 19 amino acids long
                            if i + 19 <= len(frame):
                                yield frame[i:i+19]
            else: # it's a protein, coding sequence, or other already-translated amino acid sequence, so we look at it directly
                for i in range(len(str_in)):
                    if i + 19 <= len(str_in):
                        yield str_in[i:i+19]

        # tells the interpreter that this set will not change from now on to 1) avoid bugs where we accidentally change it
        # after this point 2) optimize the data structure that is copied to every sub-process in future steps
        hazard_set = frozenset(hazard_set)
        
        # this is the function that will be run in multiple sub-processes in parallel using pool.imap_unordered() below
        def known_seqs_in_hazard_set(known_seq_fname):
            intersect_set = set()
            seq_list = []
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

        if role == 'aggregate':
            print('running as aggregator.')
            def is_partial_result(p):
                if os.path.basename(p).split('part_')[0]:
                    return False
                return os.path.basename(p).split('.')[-1] == 'json'
            indiv_paths = [p for p in indiv_paths if is_partial_result(p)]
            if not indiv_paths:
                raise RuntimeError('aggregator did not find any partial results')
            for part in indiv_paths:
                print('loading', part)
                known_hazard_seqs = set()
                with open(part) as part_f:
                    known_hazard_seqs.update(json.loads(part_f.read()))
                print('aggregated', len(known_hazard_seqs), 'known fragments from', len(indiv_paths), 'partial result files')
        else:
            print('not an aggregator; running pre-screen on known sequences in ', known_seq_path)
            print('initial size of hazard set:', len(hazard_set))
            with multiprocessing.Pool(processes=cpus) as pool:
                known_hazard_seqs = set.union(*pool.imap_unordered(known_seqs_in_hazard_set, indiv_paths))
            print('found', len(known_hazard_seqs), 'overlapping known hazard sequence fragment(s)')
            if len(known_hazard_seqs) < 100:
                print(known_hazard_seqs)
        if role == 'helper':
            print('running as helper')
            part_file_path = os.path.join(os.path.dirname(known_seq_path), 'part_' + os.path.basename(known_seq_path) + '.json')
            with open(part_file_path, 'w+') as f:
                f.write(json.dumps(list(known_hazard_seqs), indent=2))
        else:
            hazard_set = (hazard_set - known_hazard_seqs) | original_frags
            print('final size of hazard set:', len(hazard_set))
            print('writing output...')
            with open(output_path, 'w+') as f:
                f.write(json.dumps(list(hazard_set), indent=2))
    print('done')

