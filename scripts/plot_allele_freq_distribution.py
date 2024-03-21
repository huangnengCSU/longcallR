import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import sys
import gzip

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
            elif len(alt)==2 and alt[0] in ["A", "C", "G", "T"] and alt[1] in ["A", "C", "G", "T"]:
                if alt[0] == ref:
                    rna_allele_freq_dict[(chrom, pos)] = int(alt_cnt[1])/(int(alt_cnt[0])+int(alt_cnt[1]))
                elif alt[1] == ref:
                    rna_allele_freq_dict[(chrom, pos)] = int(alt_cnt[0])/(int(alt_cnt[0])+int(alt_cnt[1]))
                else:
                    ## 1/2
                    continue
    return rna_allele_freq_dict

def plot_allele_freq(rna_allele_freq_dict):
    allele_freq = list(rna_allele_freq_dict.values())
    plt.hist(allele_freq, bins=100)
    plt.xlabel('RNA allele frequency')
    plt.ylabel('Number of SNPs')
    plt.title('RNA allele frequency distribution')
    plt.savefig('allele_freq.png')

dna_snp_truth = sys.argv[1] ## DNA SNP truth set
rna_allele_freq = sys.argv[2]   ## htsbox mpileup output (htsbox pileup -f -b)
min_allele_cnt = int(sys.argv[3])   ## ignore sites with allele count less than this value
dna_hete_snps = load_dna_hete_snps(dna_snp_truth)
rna_allele_freq_dict = load_rna_allele_freq(rna_allele_freq, dna_hete_snps, min_allele_cnt)
plot_allele_freq(rna_allele_freq_dict)
