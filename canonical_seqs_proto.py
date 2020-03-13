BASES = tuple('ATCG')
from itertools import product, permutations
import time

def remap(seq, template):
    permute_map = {t:b for b, t in zip(BASES, template)}
    return tuple((permute_map[b] for b in seq))

all_pairs = [(a, b) for a in BASES for b in BASES]
pair_ids = {bp:i for i, bp in enumerate(all_pairs)}
print(pair_ids)
unordered_pair_ids = {(a, b):pair_ids[(a, b)] if BASES.index(a) > BASES.index(b) else pair_ids[(b, a)] for a, b in all_pairs}
print(unordered_pair_ids)
def order_insensitive(bps):
    return tuple(unordered_pair_ids[bp] for bp in bps)

def count_template(seq):
    t = tuple(sorted(BASES, key=lambda b:-seq.count(b)))
    print(t)
    return t

def order_template(seq):
    t = tuple(sorted(BASES, key=lambda b:seq.index(b) if b in seq else 200))
    #print(t)
    return t

    #print(canon_remap)
    #print(canon_remap)
    return canon_remap

def canon_seq(seq):
    print(seq)
    seq = remap(seq, count_template(seq))
    print(seq)
    seq = remap(seq, order_template(seq))
    print(seq)
    print('')
    halfway = len(seq)//2 # will throw out middle base if odd, fine
    return order_insensitive(zip(seq[:halfway], reversed(seq)))

print(canon_seq('AAAAAGGGG'))
print(canon_seq('GGGGAAAAA'))

def test_all(seq_len):
        for seq1 in product(BASES, repeat=n):
            for seq2 in product(BASES, repeat=n):
                is_equivalent = False
                for perm in permutations(BASES):
                    for reversal in seq1, tuple(reversed(seq1)):
                        remapped = remap(reversal, perm)
                        if seq2 == remapped:
                            is_equivalent = True
                            break
                    if is_equivalent:
                        break
                print(''.join(seq1), 'is the same as' if is_equivalent else 'is different from', ''.join(seq2))
                if is_equivalent != (canon_seq(seq1) == canon_seq(seq2)):
                    print('but we got these canon seqs:', canon_seq(seq1), canon_seq(seq2))

for n in range(3, 5):
    test_all(n)

