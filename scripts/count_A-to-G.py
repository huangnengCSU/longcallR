import sys
import gzip

hap_vcf = sys.argv[1]  ## hap.py output vcf file (decompressed)

fp_out = open("fp_sites.csv", 'w')
fn_out = open("fn_sites.csv", 'w')
non_A_to_G_tp, non_A_to_G_fp, non_A_to_G_fn = 0, 0, 0
gz_flag = False
if hap_vcf.endswith('.gz'):
    gz_flag = True
    f = gzip.open(hap_vcf, 'rb')
else:
    f = open(hap_vcf, 'r')

with f:
    tp_cnt, fp_cnt, fn_cnt = 0, 0, 0
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
            if (ref == 'A' and alt == 'G') or (ref == 'T' and alt == 'C'):
                if "TP" in truth and "SNP" in truth and "TP" in query:
                    tp_cnt += 1
                else:
                    ## An error genotype will both increase a FP and a FN
                    if "FN" in truth and "SNP" in truth:
                        fn_cnt += 1
                    if "FP" in query and "SNP" in query:
                        fp_cnt += 1
            else:
                ## non RNA-editing sites
                if "TP" in truth and "SNP" in truth and "TP" in query:
                    non_A_to_G_tp += 1
                else:
                    ## An error genotype will both increase a FP and a FN
                    if "FN" in truth and "SNP" in truth:
                        non_A_to_G_fn += 1
                        fn_out.write(line_str)
                    if "FP" in query and "SNP" in query:
                        non_A_to_G_fp += 1
                        fp_out.write(line_str)
print(f"A-to-G TP: {tp_cnt}, A-to-G FN: {fn_cnt}, A-to-G FP: {fp_cnt}")
print(f"Non A-to-G TP = {non_A_to_G_tp},\nNon A-to-G FN = {non_A_to_G_fn},\nNon A-to-G FP = {non_A_to_G_fp}")
P = non_A_to_G_tp / (non_A_to_G_tp + non_A_to_G_fp)
R = non_A_to_G_tp / (non_A_to_G_tp + non_A_to_G_fn)
F1 = 2 * P * R / (P + R)
print(f"Recall = {R},\nPrecision = {P},\nF1 = {F1}")
fp_out.close()
fn_out.close()
