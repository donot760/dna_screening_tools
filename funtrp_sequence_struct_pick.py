import os
import csv
import numpy as np

FUNTRP_DATA_DIR = os.path.join('resources', 'select_agent_hazard', 'funtrp_results')

def all_funtrp_csvs():
    return [os.path.join(FUNTRP_DATA_DIR, csv_fname)
            for csv_fname in os.listdir(FUNTRP_DATA_DIR)
            if csv_fname.split('.')[-1] == 'csv']

def csv_cols(fname, indices):
    with open(fname) as csv_file:
        reader = csv.reader(csv_file)
        csv_rows = list(reader)
    cols = []
    for col_i in indices:
        cols.append([csv_row[col_i] for csv_row in csv_rows[1:]])
    return cols

def split_into_windows(list_in, window=19):
    windows = []
    for i in range(len(list_in)):
        if i + window <= len(list_in):
            windows.append(list_in[i:i+window])
    return windows

def windows_to_floats(windows):
    return [[float(n) for n in window] for window in windows]

def add_finished_lines(lines, fragments, packed_frag_values):
    for frag in fragments:
        _, funtrp_scores = packed_frag_values[frag]
        neutrals, toggles, rheos = funtrp_scores
        lines.append(frag + ' ' + str(neutrals) + ' ' + str(toggles) + ' ' + str(rheos))

if __name__ == '__main__':
    residue_idx = 11
    neutral_idx = 6
    toggle_idx = 5
    rheo_idx = 9
    get_idxs = residue_idx, neutral_idx, toggle_idx, rheo_idx

    lines = []
    for funtrp_csv in all_funtrp_csvs():
        cols = csv_cols(funtrp_csv, get_idxs)
        windowed_cols = [split_into_windows(col) for col in cols]

        residue_windows, *float_cols = windowed_cols
        fragments = [''.join(residue_window) for residue_window in residue_windows]

        neutral_windows, toggle_windows, rheo_windows = \
                [windows_to_floats(col) for col in float_cols]
        product_scores = [np.product(np.array(toggle_window) + .3)
                for toggle_window in toggle_windows]
        packed_frag_values = {frag:(score, (neutral_window, toggle_window, rheo_window))
              for frag, score, neutral_window, toggle_window, rheo_window
              in zip(fragments, product_scores, neutral_windows, toggle_windows, rheo_windows)}
        num_frags = 10
        score_sorting_key = lambda frag: packed_frag_values[frag][0]
        sorted_fragments = list(reversed(sorted(fragments, key=score_sorting_key)))
        top_fragments = sorted_fragments[-num_frags:]
        add_finished_lines(lines, top_fragments, packed_frag_values)
    with open('resources/aa_fragment_picks.txt', 'w+') as out_f:
        for line in lines:
            out_f.write(line + '\n')

