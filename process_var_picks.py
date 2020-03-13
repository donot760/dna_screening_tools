import csv
from sequence_variants import aa_variants, m13

neutral_col = 'Select 2 neutral picks (downselect to least unusual amino acids based on BLOSUM62 diagonal)'
disruptive_col = 'Select 2 disruptive picks'

def make_ordinal(n):
    '''
    Convert an integer into its string ordinal representation

        make_ordinal(0)   => '0th'
        make_ordinal(3)   => '3rd'
        make_ordinal(122) => '122nd'
        make_ordinal(213) => '213th'
    '''
    n = int(n)
    suffix = ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)]
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    return str(n) + suffix

with open('fragments.csv') as fi, open('output_fragment_list.csv', 'w+') as fo:
    reader = csv.DictReader(fi)
    fo.write('Fragment sequence,comment')
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
        variant_map = {}
        comments = []
        for i, row in enumerate(fragment_rows):
            sub_list = [a.replace(',', '').strip() for a in row[neutral_col].split()]
            aa = row['Res']
            if sub_list:
                variant_map[(i, aa)] = sub_list
                comments.append(' and '.join(sub_list) + ' are substituted for ' + aa + ' at ' + make_ordinal(i+1) + ' position')
        comment = ', '.join(comments) + ' in ' + aa_str + ' templated after its appearance in M13'
        print(variant_map, comment)
        variants = aa_variants(aa_str, variant_map, m13)
        for seq in variants:
            fo.write(seq + ', ' + comment + '\n')
    
