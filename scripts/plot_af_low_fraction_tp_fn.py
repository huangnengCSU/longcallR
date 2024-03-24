import gzip
import sys

import matplotlib

matplotlib.use('AGG')
import matplotlib.pyplot as plt


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
                    rna_allele_freq_dict[(chrom, pos)] = [1.0, int(alt_cnt[0])]
                elif alt[0] in ["A", "C", "G", "T"] and alt[0] == ref:
                    rna_allele_freq_dict[(chrom, pos)] = [0.0, 0]
            elif len(alt) == 2 and alt[0] in ["A", "C", "G", "T"] and alt[1] in ["A", "C", "G", "T"]:
                if alt[0] == ref:
                    rna_allele_freq_dict[(chrom, pos)] = [int(alt_cnt[1]) / (int(alt_cnt[0]) + int(alt_cnt[1])), total_allele_cnt]  ## allele fraction, total allele count
                elif alt[1] == ref:
                    rna_allele_freq_dict[(chrom, pos)] = [int(alt_cnt[0]) / (int(alt_cnt[0]) + int(alt_cnt[1])), total_allele_cnt]  ## allele fraction, total allele count
                else:
                    ## 1/2
                    continue
    return rna_allele_freq_dict


def load_hap(hap_vcf):
    gz_flag = False
    if hap_vcf.endswith(".gz"):
        import gzip
        gz_flag = True
        f = gzip.open(hap_vcf, 'rb')
    else:
        f = open(hap_vcf, 'r')
    tp_cnt, fp_cnt, fn_cnt = 0, 0, 0
    tp_snps = set()
    fp_snps = set()
    fn_snps = set()
    with f:
        for line_str in f:
            if gz_flag:
                line_str = line_str.decode('utf-8')
            if line_str.startswith('#'):
                continue
            else:
                line = line_str.strip().split('\t')
                ref = line[3]
                alt = line[4]
                truth = line[9]
                query = line[10]
                if "TP" in truth and "SNP" in truth and "TP" in query:
                    tp_cnt += 1
                    tp_snps.add((line[0], int(line[1])))
                else:
                    ## An error genotype will both increase a FP and a FN
                    if "FN" in truth and "SNP" in truth:
                        fn_cnt += 1
                        fn_snps.add((line[0], int(line[1])))
                    if "FP" in query and "SNP" in query:
                        fp_cnt += 1
                        fp_snps.add((line[0], int(line[1])))
    return tp_snps, fp_snps, fn_snps


if __name__ == '__main__':
    dna_snp_truth = sys.argv[1]  ## DNA SNP truth set
    rna_allele_freq = sys.argv[2]  ## htsbox mpileup output (htsbox pileup -f -b)
    min_allele_cnt = int(sys.argv[3])  ## ignore sites with allele count less than this value
    hap_vcf = sys.argv[4]  ## hap.py output
    dna_hete_snps = load_dna_hete_snps(dna_snp_truth)
    rna_allele_freq_dict = load_rna_allele_freq(rna_allele_freq, dna_hete_snps, min_allele_cnt)
    tp_snps, fp_snps, fn_snps = load_hap(hap_vcf)

    tp_allele_fraction_list = []
    for site in tp_snps:
        if site in rna_allele_freq_dict:
            tp_allele_fraction_list.append(rna_allele_freq_dict[site][0])

    fn_allele_fraction_list = []
    for site in fn_snps:
        if site in rna_allele_freq_dict:
            fn_allele_fraction_list.append(rna_allele_freq_dict[site][0])

    plt.figure()
    plt.hist(tp_allele_fraction_list, bins=100)
    plt.xlabel('Allele fraction')
    plt.ylabel('Number of SNPs')
    plt.title('Allele fraction distribution of TP low allele fraction SNPs')
    plt.savefig('tp_low_frac.png')

    plt.figure()
    plt.hist(fn_allele_fraction_list, bins=100)
    plt.xlabel('Allele fraction')
    plt.ylabel('Number of SNPs')
    plt.title('Allele fraction distribution of FN low allele fraction SNPs')
    plt.savefig('fn_low_frac.png')


    plt.figure()
    plt.hist(tp_allele_fraction_list, bins=100, alpha = 0.4, color = 'red')
    plt.hist(fn_allele_fraction_list, bins=100, alpha = 0.4, color = 'blue')
    plt.legend(['tp snps', 'fn snps'])
    plt.xlabel('Allele fraction')
    plt.ylabel('Number of SNPs')
    plt.title('Allele fraction distribution of TP and FN low allele fraction SNPs')
