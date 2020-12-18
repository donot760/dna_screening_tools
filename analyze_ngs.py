from Bio import SeqIO
from Bio.Seq import Seq
import pdb
import csv
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from random import sample, random

from ecoli_seq_variants import single_variants
from funtrp_blosum_predict import fragment_replace_mtx
from metro_hastings_variants import int_encode_variant
from fastq_counting import fastq_fragment_counts

data_dir = os.path.join('resources', 'NGS', '6_27_DSS_libs_round1')
#file_names = ["ESP7xLib14-preselection-miniprep.fastq"]
file_names = ["ESP7xLib14-postselection-phage.fastq"]
#file_names = ["ESP7xLib14-preselection-miniprep.fastq", "ESP7xLib14-postselection-phage.fastq"]
fragment = sys.argv[sys.argv.index('--fragment')+1] if '--fragment' in sys.argv else 'YANYEGCLWNATGVVVCTG'

## Read the fastq files. Count the total number of sequences
all_data = {}
num_entries = 0
def load_all_data():
    global all_data, num_entries
    for f in file_names:
        all_data[f] = []
        with open(os.path.join(data_dir, f)) as f_handle:
            for record in SeqIO.parse(f_handle, "fastq"):
                num_entries += 1
                all_data[f].append(record.seq)
        print("Total reads", f, num_entries)

def pad_codons(seq):
    return seq + (-len(seq) % 3)*'N'

def translate(seq, mem={}):
    seq = seq.lower()
    if seq in mem:
       return mem[seq]
    padded = pad_codons(seq)
    try:
        translation = str(padded.translate())
    except (AttributeError, TypeError):
        padded = Seq(padded)
        translation = str(padded.translate())
    mem[seq] = translation
    return translation

def count_stats(seq):
    count_exact = 0
    count_trans = 0
    seq = seq.lower()
    trans = translate(seq)
    for fname in file_names:
        for single_read in all_data[fname]:
            rev_comp = single_read.reverse_complement()
            #print(seq, 'in', single_read, 'or', rev_comp, '?')
            if seq in str(single_read).lower() or seq in str(rev_comp).lower():
                count_exact += 1
            for frame in range(3):
                read_trans = translate(single_read[frame:])
                read_comp_trans = translate(rev_comp[frame:])
                if trans in read_trans or trans in read_comp_trans:
                    count_trans += 1
    assert count_trans >= count_exact
    return count_exact, count_trans

wt_sequence = 'tacgctaactatgagggctgtctgtggaatgctacaggcgttgtagtttgtactggt'.upper()
wt_trans = translate(wt_sequence)
#wt_sequence = 'tacgctaactatgagggc'.upper()
# wt translation: YANYEGCLWNATGVVVCTG

single_vars = None
def load_single_vars():
    global single_vars
    single_vars = single_variants(wt_trans, orig_dna_seq=wt_sequence)


all_vars = all_transs = original_fragments = None
def load_all_vars_and_transs():
    global all_vars, all_transs, original_fragments
    with open('/Users/danagretton/Dropbox (MIT)/Sculpting Evolution/DNA Sequence Validation/Screening M13/Quotes Round 2/extracted_submitted_aa_variants.json') as f:
        all_vars = set(json.loads(f.read())) - {'', '\n', '\r\n'}
    all_transs = set((translate(var) for var in all_vars))
    with open('/Users/danagretton/Dropbox (MIT)/Sculpting Evolution/DNA Sequence Validation/Screening M13/Quotes Round 2/original_fragment_map.json') as f:
        original_fragments = json.loads(f.read())
    all_transs = set((v for v in all_transs if v == fragment or original_fragments.get(v, None) == fragment))
    print('number of transs:', len(all_transs))

def windows(seq):
    for i in range(len(seq)):
        window = seq[i:i+57]
        if len(window) != 57:
            return
        yield window.lower()

# count instances of this sequence
trans_counts = {}
def count_translations():
    count_wt, count_wt_trans = count_stats(wt_sequence)
    print("WT exact sequence count", count_wt)
    print("WT synonymous sequence count", count_wt_trans)
    for variant in all_vars:#[:40]:
        count_vari, count_vari_trans = count_stats(variant)
        trans_counts[translate(variant)] = count_vari_trans
        if count_vari_trans > 0:
            print(variant)
            print(translate(variant))
            print("Variant exact sequence count", count_vari)
            print("Variant synonymous sequence count", count_vari_trans)

count_is_wt = count_is_exact = count_is_aa_match = count_is_other = 0
covered_vars = set()
covered_transs = set()
def add_single_read(single_read):
    global count_is_wt, count_is_exact, count_is_aa_match, count_is_other
    rev_comp = single_read.reverse_complement()
    #print(seq, 'in', single_read, 'or', rev_comp, '?')
    seq = wt_sequence
    seq = seq.lower()
    if seq in str(single_read).lower() or seq in str(rev_comp).lower():
        print('found wild type')
        count_is_wt += 1
        return
    for window in windows(str(single_read)):
        if window in all_vars:
            covered_vars.add(window)
            covered_transs.add(window)
            count_is_exact += 1
            return
        var_translation = translate(window)
        if var_translation in all_transs:
            covered_transs.add(window)
            count_is_aa_match += 1
            return
    for window in windows(str(rev_comp)):
        if window in all_vars:
            covered_vars.add(window)
            covered_transs.add(window)
            count_is_exact += 1
            return
        var_translation = translate(window)
        if var_translation in all_transs:
            covered_transs.add(window)
            count_is_aa_match += 1
            return
    count_is_other += 1

replace_mtx = fragment_replace_mtx(fragment)
log_replace_mtx = np.log(replace_mtx)
def variant_log_prob(variant):
    return sum((log_replace_mtx[i,j] for j, i in enumerate(int_encode_variant(variant))))

def percent(c, denom=num_entries):
    return str(c/denom*100) + '%'

cache_dir = os.path.join('resources', 'cache')
def get_counts(fnames):
    combined_name = os.path.join(cache_dir, '--'.join(sorted((
            [os.path.basename(fname).split('.')[0] for fname in fnames]
            ))) + '.cache.json')
    try:
        with open(combined_name) as cache_f:
            return json.loads(cache_f.read())
    except FileNotFoundError:
        counts = fastq_fragment_counts(fnames, leading_seq='tttagatcgt',
                trailing_seq=None)#'gacgaaactc')
        with open(combined_name, 'w+') as cache_f:
            cache_f.write(json.dumps(counts, indent=4))
        return counts

if __name__ == '__main__':
    print('Fragment:' + fragment)
    #load_all_data()
    #count_translations()
    load_all_vars_and_transs()
    import matplotlib.pyplot as plt

    #for single_read in all_data["ESP7xLib14-postselection-phage.fastq"]:
    #    add_single_read(single_read)

    ngs_dir = os.path.join('resources', 'NGS', '20200810 Genewiz Results', 'Flattened and Expanded NGS Results')
    pre_counts = get_counts([os.path.join(ngs_dir, '30-395819904__00_fastq__22-7-24_R1_001.fastq'), os.path.join(ngs_dir, '30-395819904__00_fastq__22-7-24_R2_001.fastq')])
    post_counts = get_counts(['/Users/danagretton/Dropbox (MIT)/Sculpting Evolution/Dev/dna_screening_tools/resources/NGS/6_27_DSS_libs_round1/ESP7xLib14-postselection-phage.fastq'])
    xs = []
    ys = []
    #for trans in set(pre_counts) | set(post_counts):
    for trans in all_transs:
        #if trans in pre_counts and trans in post_counts:
        #    xs.append(pre_counts[trans])
        #    ys.append(post_counts[trans])
        x = pre_counts.get(trans, 0)
        y = post_counts.get(trans, 0)
        if x > 140000 and 100 < y < 6000:
            print(trans)
        if x > 0 and y > 0:
            print(trans)
        if x + y > 500:
            xs.append(x)
            ys.append(y)# + random()*.15)
    import matplotlib.pyplot as plt
    #plt.scatter(xs, ys)
    #plt.show()

    plot_variants = []
    for i, var in enumerate(post_counts):
        if '*' in var:
            continue
        if var not in pre_counts:
            continue
        plot_variants.append(var)
    pre_norm = sum(pre_counts.values())
    post_norm = sum(post_counts.values())
    xs = [np.log(post_counts[var]/post_norm/(pre_counts.get(var, 1)/pre_norm)) for var in plot_variants]
    ys = [np.exp(variant_log_prob(var)) for var in plot_variants]
    plt.scatter(xs, ys)
    for var, x, y in zip(plot_variants, xs, ys):
        label = var + ": {:.5f}".format(x) + (' (WILD TYPE)' if var == fragment else '')

        # this method is called for each point
        plt.annotate(label, # this is the text
                     (x,y), # this is the point to label
                     textcoords="offset points", # how to position the text
                     xytext=(0,10), # distance from text to points (x,y)
                     ha='center') # horizontal alignment can be left, right or center
    plt.show()

    import pdb; pdb.set_trace()

    #print(percent(count_is_wt), percent(count_is_exact), percent(count_is_aa_match), percent(count_is_other))
    #print('total:', percent(count_is_wt + count_is_exact + count_is_aa_match + count_is_other))
    #print('exact library sequence coverage:', percent(len(covered_vars), len(all_vars)))
    #print('aa library coverage:', percent(len(covered_transs), len(all_transs)))

