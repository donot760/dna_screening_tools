from Bio import SeqIO
from Bio.Seq import Seq
import pdb
import csv
import os
import json

from ecoli_seq_variants import single_variants

data_dir = '6_27_DSS_libs_round1'
#file_names = ["ESP7xLib14-preselection-miniprep.fastq"]
file_names = ["ESP7xLib14-postselection-phage.fastq"]
#file_names = ["ESP7xLib14-preselection-miniprep.fastq", "ESP7xLib14-postselection-phage.fastq"]

## Read the fastq files. Count the total number of sequences
all_data = {}
for f in file_names:
    all_data[f] = []
    num_entries = 0
    with open(os.path.join(data_dir, f)) as handle:
        for record in SeqIO.parse(handle, "fastq"):
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

def count_stats(fname, seq):
    count_exact = 0
    count_trans = 0
    seq = seq.lower()
    trans = translate(seq)
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

single_vars = single_variants(wt_trans, orig_dna_seq=wt_sequence)
with open('/Users/danagretton/Dropbox (MIT)/Sculpting Evolution/DNA Sequence Validation/Screening M13/Quotes Round 2/extracted_submitted_aa_variants.json') as f:
    all_vars = set(json.loads(f.read()))
all_transs = set((translate(var) for var in all_vars))

def windows(seq):
    for i in range(len(seq)):
        window = seq[i:i+57]
        if len(window) != 57:
            return
        yield window.lower()

# count instances of this sequence
for fname in []:#file_names: #TODO: this is disabled
    print(fname)
    count_wt, count_wt_trans = count_stats(fname, wt_sequence)
    print("WT exact sequence count", fname, count_wt)
    print("WT synonymous sequence count", fname, count_wt_trans)
    for variant in single_vars:#[:40]:
        count_vari, count_vari_trans = count_stats(fname, variant)
        if count_vari_trans > 0:
            print(variant)
            print(translate(variant))
            print("Variant exact sequence count", fname, count_vari)
            print("Variant synonymous sequence count", fname, count_vari_trans)

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
        if translate(window) in all_transs:
            covered_transs.add(window)
            count_is_aa_match += 1
            return
    for window in windows(str(rev_comp)):
        if window in all_vars:
            covered_vars.add(window)
            covered_transs.add(window)
            count_is_exact += 1
            return
        if translate(window) in all_transs:
            covered_transs.add(window)
            count_is_aa_match += 1
            return
    count_is_other += 1

for single_read in all_data["ESP7xLib14-postselection-phage.fastq"]:
    add_single_read(single_read)

def percent(c, denom=num_entries):
    return str(c/denom*100) + '%'

print(percent(count_is_wt), percent(count_is_exact), percent(count_is_aa_match), percent(count_is_other))
print('total:', percent(count_is_wt + count_is_exact + count_is_aa_match + count_is_other))
print('exact library sequence coverage:', percent(len(covered_vars), len(all_vars)))
print('aa library coverage:', percent(len(covered_transs), len(all_transs)))

