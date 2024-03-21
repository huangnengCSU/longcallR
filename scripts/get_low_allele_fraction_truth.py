import gzip
import sys


def load_dna_hete_snps(dna_snp_truth):
    dna_hete_snps = set()
    gz_flag = False
    if dna_snp_truth.endswith(".gz"):
        gz_flag = True
        f = gzip.open(dna_snp_truth, 'rb')
    else:
        f = open(dna_snp_truth, 'r')
    for line in f:
        if gz_flag:
            line = line.decode('utf-8')
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        gt = fields[-1].split(":")[0]
        if gt == "0/1" or gt == "1/0":
            dna_hete_snps.add((chrom, pos))
    f.close()
    return dna_hete_snps


def load_rna_allele_freq(rna_allele_freq, dna_hete_snps, min_allele_cnt):
    rna_allele_freq_dict = {}
    with open(rna_allele_freq, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            if (chrom, pos) not in dna_hete_snps:
                continue
            ref = fields[2]
            alt = fields[3].split(",")
            alt_cnt = fields[4].split(":")[1].split(",")
            total_allele_cnt = sum([int(x) for x in alt_cnt])
            if total_allele_cnt < min_allele_cnt:
                continue
            if len(alt) == 1:
                if alt[0] in ["A", "C", "G", "T"] and alt[0] != ref:
                    rna_allele_freq_dict[(chrom, pos)] = 1.0
                elif alt[0] in ["A", "C", "G", "T"] and alt[0] == ref:
                    rna_allele_freq_dict[(chrom, pos)] = 0.0
            elif len(alt) == 2 and alt[0] in ["A", "C", "G", "T"] and alt[1] in ["A", "C", "G", "T"]:
                if alt[0] == ref:
                    rna_allele_freq_dict[(chrom, pos)] = int(alt_cnt[1]) / (int(alt_cnt[0]) + int(alt_cnt[1]))
                elif alt[1] == ref:
                    rna_allele_freq_dict[(chrom, pos)] = int(alt_cnt[0]) / (int(alt_cnt[0]) + int(alt_cnt[1]))
                else:
                    ## 1/2
                    continue
    return rna_allele_freq_dict


def output_low_allele_fraction_truth(dna_snp_truth, rna_allele_freq_dict, min_allele_fraction, max_allele_fraction,
                                     output_file):
    low_allele_fraction_sites = set()
    for k, v in rna_allele_freq_dict.items():
        if v < max_allele_fraction and v >= min_allele_fraction:
            low_allele_fraction_sites.add(k)

    fout = open(output_file, 'w')
    gz_flag = False
    if dna_snp_truth.endswith(".gz"):
        gz_flag = True
        f = gzip.open(dna_snp_truth, 'rb')
    else:
        f = open(dna_snp_truth, 'r')
    for line in f:
        if gz_flag:
            line = line.decode('utf-8')
        if line.startswith("#"):
            fout.write(line)
            continue
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        if (chrom, pos) in low_allele_fraction_sites:
            fout.write(line)
    f.close()
    fout.close()


if __name__ == "__main__":
    dna_snp_truth = sys.argv[1]  ## DNA SNP truth set
    rna_allele_freq = sys.argv[2]  ## htsbox mpileup output (htsbox pileup -f -b)
    min_allele_cnt = int(sys.argv[3])  ## ignore sites with allele count less than this value
    min_allele_fraction = float(sys.argv[4])
    max_allele_fraction = float(sys.argv[5])
    output_file = sys.argv[6]
    dna_hete_snps = load_dna_hete_snps(dna_snp_truth)
    rna_allele_freq_dict = load_rna_allele_freq(rna_allele_freq, dna_hete_snps, min_allele_cnt)
    output_low_allele_fraction_truth(dna_snp_truth, rna_allele_freq_dict, min_allele_fraction, max_allele_fraction,
                                     output_file)
