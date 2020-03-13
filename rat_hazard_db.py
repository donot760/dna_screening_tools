from Bio.Seq import Seq
from hashlib import sha256
from hazard_db import split_into_windows, hazard_db_hash, canon_seq
from hazard_db import SEQ_WINDOW, AMINO_WINDOW, M13
import numpy as np
import json

aas = 'ARNDCQEGHILKMFPSTWYV'

aa_replace_mtx = np.array([
[0,0,0,0,0,0,1,3,0,0,0,0,0,0,2,4,3,0,0,1],
[0,0,0,0,0,4,0,0,3,0,0,6,0,0,0,0,0,2,0,0],
[0,0,0,6,0,0,2,0,4,0,0,3,0,0,0,3,2,0,0,0],
[0,0,6,0,0,1,7,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,4,1,1,0,0,6,0,6,0,0,3,1,0,0,0,0,0,0,0],
[1,0,2,7,0,6,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
[3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
[0,3,4,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,2,0],
[0,0,0,0,0,0,0,0,0,0,4,0,4,1,0,0,1,0,0,8],
[0,0,0,0,0,0,0,0,0,4,0,0,7,2,0,0,0,0,0,2],
[0,6,3,0,0,3,1,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,4,7,1,0,0,0,0,0,0,0,3],
[0,0,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,1,9,0],
[2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
[4,0,3,0,0,0,0,1,0,0,0,0,0,0,1,0,5,0,0,0],
[3,0,2,0,0,0,0,0,0,1,0,0,0,0,0,5,0,0,0,1],
[0,2,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,3,0],
[0,0,0,0,0,0,0,0,2,0,0,0,0,9,0,0,0,3,0,0],
[1,0,0,0,0,0,0,0,0,8,2,0,3,0,0,0,1,0,0,0]
])

def select_windows(window_len, seq, dist, frac):
    # returns a list of string windows of length window_len from string seq sampled
    # according to a distribution dist, list of len(seq) - window_len + 1 floats
    # summing to 1, until frac of the sequence is covered
    all_windows = split_into_windows(seq, window_len)
    dist = dist[window_len//2:(-window_len)//2+1] # use centers of windows
    dist = np.array(dist)/sum(dist)
    windows = set()
    while len(windows) < len(all_windows)*frac:
        windows.add(np.random.choice(all_windows, p=dist))
    return windows
    all_windows = split_into_windows(seq, window_len)
    dist = dist[window_len//2:(-window_len)//2+1] # use centers of windows
    dist = np.array(dist)/sum(dist)
    windows = set()
    while len(windows) < len(all_windows)*frac:
        windows.add(np.random.choice(all_windows, p=dist))
    return windows

def expand_window(aa_window, depth, num_changes=1):
    # return every substitution scoring at least depth (int 0-10) of a string of
    # amino acid abbreviation letters
    if num_changes != 1:
        raise NotImplementedError
    variations = [aa_window]
    for i, aa in enumerate(aa_window):
        for aa2 in aas:
            try:
                j_start, j_replace = aas.index(aa), aas.index(aa2)
            except ValueError:
                continue
            if 9 - aa_replace_mtx[j_start, j_replace] <= depth:
                variations.append(aa_window[:i] + aa2 + aa_window[i+1:])
    return variations

class RAT_db:
    def __init__(self):
        self.db = set()

    def save(self, path):
        with open(path, 'w+') as f:
            f.write(json.dumps(list(self.db)))

    def load(self, path):
        with open(path) as f:
            for entry in json.loads(f.read()):
                self.db.add(entry)

    def add_file(self, path, conserve_dist, depth, cover_frac):
        with open(path) as f:
            seq = f.read().strip()
        self.add(seq, conserve_dist, depth, cover_frac)

    def screen_file(self, path):
        with open(path) as f:
            seq = f.read().strip()
        return self.screen(seq)

    def add(self, seq, conserve_dist, depth, cover_frac):
        for seq_window in split_into_windows(seq, SEQ_WINDOW):
            self.db.add(hazard_db_hash(canon_seq(seq_window)))
        for aa_window in select_windows(AMINO_WINDOW*3, seq, conserve_dist, cover_frac):
            translation = str(Seq(aa_window).translate())
            for variation in expand_window(translation, depth):
                self.db.add(hazard_db_hash(variation))

    def screen(self, seq):
        for window in split_into_windows(seq, SEQ_WINDOW):
            if hazard_db_hash(canon_seq(window)) in self.db:
                return False
        for seq_dir in seq, str(Seq(seq).reverse_complement()):
            for aa_window in split_into_windows(seq_dir, AMINO_WINDOW*3):
                translation = str(Seq(aa_window).translate())
                if hazard_db_hash(translation) in self.db:
                    return False
        return True

if __name__ == '__main__':
    window, *_ = split_into_windows(M13, SEQ_WINDOW)
    print(expand_window(str(Seq(window).translate()), 5))
    rdb = RAT_db()
    cover_frac = .4
    rdb.add(M13, [np.random.random() for b in M13], 2, cover_frac)
    for window in split_into_windows(M13, SEQ_WINDOW):
        if hazard_db_hash(canon_seq(window)) not in rdb.db:
            raise RuntimeError
        print(window, 'checks out')
    num_present = 0
    aa_windows = split_into_windows(M13, AMINO_WINDOW*3)
    num_expected = len(aa_windows) * cover_frac
    for aa_window in aa_windows:
        translation = str(Seq(aa_window).translate())
        if hazard_db_hash(translation) in rdb.db:
            print(translation, 'checks out')
            num_present += 1
        if num_present > num_expected:
            break
    else:
        raise RuntimeError
    assert rdb.screen(''.join(np.random.choice(list('atcg'), 100000)))
    assert not rdb.screen(M13)
    print('screened no false positives ok')
    rdb.save('m13_f_40_d_2.json')
    rdb = RAT_db()
    assert rdb.screen(M13)
    rdb.load('m13_f_40_d_2.json')
    assert not rdb.screen(M13)
    print('saved and loaded ok')
    rdb = RAT_db()
    assert rdb.screen(M13)
    assert rdb.screen_file('m13.seq')
    rdb.add_file('m13.seq', [np.random.random() for b in M13], 2, cover_frac)
    assert not rdb.screen(M13)
    assert not rdb.screen_file('m13.seq')
    print('added sequence from file ok')
    print('Everything checks out!')

