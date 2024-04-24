import matplotlib

matplotlib.use('AGG')
import matplotlib.pyplot as plt
import sys


def load_rna_edits(rna_edit_file):
    rna_edits = set()
    with open(rna_edit_file, 'r') as f:
        f.readline()
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            rna_edits.add((chrom, pos))
    return rna_edits


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
                ## homogyzous
                if alt[0] == "G" and ref == "A" and int(alt_cnt[0]) >= 2:
                    rna_allele_freq_dict[(chrom, pos)] = 1.0
                elif alt[0] == "C" and ref == "T" and int(alt_cnt[0]) >= 2:
                    rna_allele_freq_dict[(chrom, pos)] = 1.0
            elif len(alt) == 2:
                if alt[0] == "A" and alt[1] == "G" and int(alt_cnt[1]) >= 2:
                    rna_allele_freq_dict[(chrom, pos)] = int(alt_cnt[1]) / (int(alt_cnt[0]) + int(alt_cnt[1]))
                elif alt[0] == "T" and alt[1] == "C" and int(alt_cnt[1]) >= 2:
                    rna_allele_freq_dict[(chrom, pos)] = int(alt_cnt[1]) / (int(alt_cnt[0]) + int(alt_cnt[1]))
                elif alt[0] == "G" and alt[1] == "A" and int(alt_cnt[0]) >= 2:
                    rna_allele_freq_dict[(chrom, pos)] = int(alt_cnt[0]) / (int(alt_cnt[0]) + int(alt_cnt[1]))
                elif alt[0] == "C" and alt[1] == "T" and int(alt_cnt[0]) >= 2:
                    rna_allele_freq_dict[(chrom, pos)] = int(alt_cnt[0]) / (int(alt_cnt[0]) + int(alt_cnt[1]))
    return rna_allele_freq_dict


def plot_allele_freq(rna_allele_freq_dict):
    allele_freq = list(rna_allele_freq_dict.values())
    plt.hist(allele_freq, bins=10)
    plt.xlabel('RNA editing allele frequency')
    plt.ylabel('Number of SNPs')
    plt.title('RNA editing allele frequency distribution')
    plt.savefig('rna_endt_allele_freq.png')


rna_edit_dataset_file = sys.argv[1]  ## RNA editing truth file
rna_allele_freq = sys.argv[2]  ## htsbox mpileup output (htsbox pileup -f -b)
min_allele_cnt = int(sys.argv[3])  ## ignore sites with allele count less than this value
dna_hete_snps = load_rna_edits(rna_edit_dataset_file)
rna_allele_freq_dict = load_rna_allele_freq(rna_allele_freq, dna_hete_snps, min_allele_cnt)
plot_allele_freq(rna_allele_freq_dict)
