from itertools import product, permutations
import time

BASES = tuple(range(4))
def str_seq(int_seq):
    return ''.join(('atcg'[i] for i in int_seq))

def to_int_seq(seq):
    return tuple((BASES['atcg'.index(b)] for b in seq.lower()))

known_maps = {}
def remap(int_seq, template):
    if template in known_maps:
        permute_map = known_maps[template]
    else:
        permute_map = {t:b for b, t in zip(BASES, template)}
        known_maps[template] = permute_map
    return tuple((permute_map[b] for b in int_seq))

def order_template(int_seq):
    return tuple(sorted(BASES, key=lambda b:int_seq.index(b) if b in int_seq else 2000))

def canon_seq(seq):
    int_seq = to_int_seq(seq)
    rev = tuple(reversed(int_seq))
    int_seq = remap(int_seq, order_template(int_seq))
    if rev < int_seq:
        int_seq = rev
    return str_seq(int_seq)

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
            if is_equivalent != (canon_seq(seq1) == canon_seq(seq2)):
                print(''.join(seq1), 'is the same as' if is_equivalent else 'is different from', ''.join(seq2))
                print('but we got these canon seqs:', canon_seq(seq1), canon_seq(seq2))
                return False
    return True

if __name__ == '__main__':
    for n in range(2, 5):
        if test_all(n):
            print('All', n, 'check out.')
        else:
            break

