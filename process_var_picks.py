import csv
from ecoli_seq_variants import aa_variants, single_variants, m13, unique, all_aminos, CommentedList
from itertools import zip_longest

neutral_col = 'Select 4 neutral picks (downselect to least unusual amino acids based on BLOSUM62 diagonal)'
disruptive_col = 'Select 2 disruptive picks'

with open('fragment_variant_picks.csv') as fi, open('aa_fragment_variants.csv', 'w+') as fo:
    reader = csv.DictReader(fi)
    fo.write('Comment,Fragment sequence\n')
    while True:
        fragment_rows = []
        aa_str = ''
        for row in reader: # max fragment len 100
            if not row['Position']: # first column empty indicates break between fragments
                break
            if not row['Position']: # eat up extra empty rows
                continue
            fragment_rows.append(row)
            aa_str += row['Res']
        if not aa_str:
            break # outer loop
        print(aa_str)
        variants = CommentedList()
        for select_col in neutral_col, disruptive_col:
            variant_map = {}
            annotations = []
            for i, row in enumerate(fragment_rows):
                sub_list = [a for a in row[select_col].replace(',', ' ').split()]
                aa = row['Res']
                annotations.append(row['Position'])
                if sub_list:
                    variant_map[(i, aa)] = sub_list
            variants.extend(aa_variants(aa_str, variant_map, m13, annotations=annotations, padding=26))
        one_variants = single_variants(aa_str, m13, annotations=['single variants of ' + a for a in annotations], padding=26)
        variants.extend(one_variants)
        if True: # whether to include all pairs of neutral variants
            for idx, aa in zip(range(len(aa_str)), aa_str):
                row = fragment_rows[idx]
                if row['PICKS'].strip() != '*':
                    continue
                for fwd_idx, fwd_aa in zip(range(len(aa_str)), aa_str):
                    if fwd_idx <= idx:
                        continue # avoid inefficient double counting
                    fwd_row = fragment_rows[fwd_idx]
                    if fwd_row['PICKS'].strip() != '*':
                        continue # from this point, only have pairs of neutral picks
                    print(aa_str)
                    print(aa, fwd_aa)
                    variant_map = {(idx, aa):all_aminos, (fwd_idx, fwd_aa):all_aminos}
                    variants.extend(aa_variants(aa_str, variant_map, m13, annotations=annotations, padding=26))
        variants = unique(variants)
        for seq, comment in zip_longest(variants, variants.comments):
            fo.write('"' + comment + '",' + seq + '\n')

