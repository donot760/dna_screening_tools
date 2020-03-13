from itertools import product
from Bio.Seq import Seq

with open('m13.seq') as f:
    m13 = f.read()

ecoli_codons = {
'A': 'gcg',
'R': 'cgc',
'N': 'aac',
'D': 'gat',
'C': 'tgc',
'E': 'gaa',
'Q': 'cag',
'G': 'ggc',
'H': 'cat',
'I': 'att',
'L': 'ctg',
'K': 'aaa',
'M': 'atg',
'F': 'ttt',
'P': 'ccg',
'S': 'agc',
'T': 'acc',
'W': 'tgg',
'Y': 'tat',
'V': 'gtg'
}

all_aminos = list(ecoli_codons.keys())

def unique(coll):
    return list(set(coll))

def aa_variants(aa_seq_in, replace_map, source_dna_seq=None):
    """
    E. coli-codon-optimized dna variants of an amino acid sequence.
    aa_seq_in: string, single-letter amino acid codes
    replace_map: dictionary mapping tuples to lists.
        - Tuples are (index in aa_seq_in, single-letter aa in aa_seq_in) for validation.
        - Lists are single-letter aa codes of replacement options.
    [source_dna_seq]: An optional dna string of arbitrary length.
        - If given, the first instance of aa_seq_in in its translation will be used to
          fill in codons for unchanged aas. All 3 frames are checked in sequence.

    returns: list of variants as [atcg]* strings
        - All replacement combinations are included, including the original aa sequence.
    """
    if source_dna_seq is None:
        orig_dna_seq = ''.join((ecoli_codons[aa] for aa in aa_seq_in))
    else:
        for offset in range(3):
            frame = source_dna_seq[offset:]
            frame = frame[:len(frame)//3*3] # trim to satisfy BioPython
            try:
                pattern_index = str(Seq(frame).translate()).index(aa_seq_in)*3
            except ValueError:
                print('Not found in this frame!!')
                continue
            orig_dna_seq = frame[pattern_index:pattern_index+len(aa_seq_in)*3]
            break
        else:
            raise ValueError('source DNA given, but aa sequence ' + aa_seq_in + ' does not appear in any frame of its translation.')
    length = len(aa_seq_in)
    assert len(orig_dna_seq)/3 == length
    for replace_idx, replace_aa in replace_map:
        # make sure identity replace is in replace map
        replace_map[(replace_idx, replace_aa)] = set(list(replace_map[(replace_idx, replace_aa)]) + [replace_aa])
    all_replaces = [[(idx_and_aa, replace_aa) for replace_aa in replace_map[idx_and_aa]] for idx_and_aa in replace_map]
    variants = []
    for replace_combo in product(*all_replaces):
        new_aa_seq = aa_seq_in
        for (replace_idx, old_aa), new_aa in replace_combo:
            if new_aa_seq[replace_idx] != old_aa:
                raise ValueError
            new_aa_seq = new_aa_seq[:replace_idx] + new_aa + new_aa_seq[replace_idx+1:]
        seq_out = ''
        for old_aa, new_aa, dna_idx in zip(aa_seq_in, new_aa_seq, range(0, len(orig_dna_seq), 3)):
            if new_aa == old_aa:
                seq_out += orig_dna_seq[dna_idx:dna_idx+3]
            else:
                seq_out += ecoli_codons[new_aa]

        variants.append(seq_out)
    return unique(variants)

def single_variants(aa_seq_in, orig_dna_seq=None):
    variants = []
    for idx_and_aa in zip(range(len(aa_seq_in)), aa_seq_in):
        replacements = {idx_and_aa:all_aminos}
        variants += aa_variants(aa_seq_in, replacements, orig_dna_seq)
    return unique(variants)

if __name__ == '__main__':
    aa_fragment = 'NSPLMNNFRQYLPSLPQSV'
    neutral_replacements = {
        (0, 'N'): ['T', 'S'],
        (1, 'S'): ['A', 'T'],
        (2, 'P'): ['A', 'T'],
        (5, 'N'): ['S', 'D'],
        (10, 'Y'): ['H', 'F'],
        (13, 'S'): ['A', 'T']
    }
    variants = sorted(unique(single_variants(aa_fragment, m13) + aa_variants(aa_fragment, neutral_replacements, m13)))
    print('Genrated', len(variants), 'variants')
    with open('variants.txt', 'w+') as f:
        for variant in variants:
            f.write(variant + '\n')

