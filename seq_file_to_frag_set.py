import sys, json
if len(sys.argv) < 3:
    print('must supply input and output file names')
    exit()
allowed_chars, window = ('ATCG', 42) if '--dna' in sys.argv else ('ARNDBCEQZGHILKMFPSTWYV', 19)
with open(sys.argv[-2]) as f:
    a = f.read().replace('\n', '').replace(' ', '')
out_set = set()
out_list = []
for i, j in zip(range(len(a)), range(window, len(a) + 1)):
    fragment = a[i:j]
    if fragment not in out_set and all((c in allowed_chars for c in fragment)):
        out_set.add(fragment)
        out_list.append(fragment)
with open(sys.argv[-1], 'w+') as f:
    f.write(json.dumps(out_list).replace('",', '",\n'))
