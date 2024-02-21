import argparse
from collections import defaultdict


def read_isoform_file(isoform_file):
    read_isoform = defaultdict(list)
    with open(isoform_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            read_isoform[line[1]].append(line[0])
    return read_isoform


def read_haplotype_file(haplotype_file):
    read_haplotype = {}
    with open(haplotype_file, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            read_haplotype[line[0]] = line[1]
    return read_haplotype


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("-i", "--isoform", type=str, required=True, help="read isoform file")
    args.add_argument("-a", "--haplotype", type=str, required=True, help="read haplotype file")
    args.add_argument("-s", "--min_support", type=int, default=8,
                      help="minimum support of reads from one haplotype for haplotype-specific isoform")
    args = args.parse_args()
    read_isoform = read_isoform_file(args.isoform)
    read_haplotype = read_haplotype_file(args.haplotype)
    for isoform_id, read_ids in read_isoform.items():
        hap1_cnt = 0
        hap2_cnt = 0
        hap1_reads = []
        hap2_reads = []
        for read_id in read_ids:
            if read_id in read_haplotype:
                if read_haplotype[read_id] == "1":
                    hap1_cnt += 1
                    hap1_reads.append(read_id)
                elif read_haplotype[read_id] == "2":
                    hap2_cnt += 1
                    hap2_reads.append(read_id)
        if hap1_cnt + hap2_cnt >= args.min_support and hap1_cnt * hap2_cnt == 0:
            print(isoform_id, hap1_cnt, hap2_cnt)
            print(hap1_reads)
            print(hap2_reads)
