import gzip
import sys


def load_rna_edit_calls(rna_edit_calls_file):
    rna_edit_calls = set()
    gz_tag = False
    if rna_edit_calls_file.endswith(".gz"):
        gz_tag = True
        f = gzip.open(rna_edit_calls_file, 'rb')
    else:
        f = open(rna_edit_calls_file, 'r')
    for line in f:
        if gz_tag:
            line = line.decode('utf-8')
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        rna_edit_calls.add((chrom, pos))
    f.close()
    return rna_edit_calls


def load_truth(truth_file):
    truth = set()
    with open(truth_file, 'r') as f:
        f.readline()    ## skip header
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            truth.add((chrom, pos))
    return truth


if __name__ == '__main__':
    rna_edit_calls_file = sys.argv[1]
    truth_file = sys.argv[2]
    rna_edit_calls = load_rna_edit_calls(rna_edit_calls_file)
    truth = load_truth(truth_file)
    tp_fp = len(rna_edit_calls)
    tp = 0
    for call in rna_edit_calls:
        if call in truth:
            tp += 1
    precision = tp / tp_fp
    print("tp: %d, fp: %d, precision: %.4f" % (tp, tp_fp - tp, precision))
