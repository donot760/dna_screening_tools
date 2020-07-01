import random

quality_encoding = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
def fastq(label, seq, qualities):
    if len(seq) != len(qualities):
        raise ValueError
    return '@' + label + '\n' + seq + '\n+\n' + ''.join((quality_encoding[q] for q in qualities))

def random_fastq(label, seq):
    return fastq(label, seq, [random.randint(0, len(quality_encodings)) for _ in seq])

