import sys
import gzip

hap_vcf = sys.argv[1]  ## hap.py output vcf file (decompressed)

fp_out = open("fp_sites.csv", 'w')
fn_out = open("fn_sites.csv", 'w')
non_A_to_I_tp, non_A_to_I_fp, non_A_to_I_fn = 0, 0, 0
non_A_to_I_tp_snps, non_A_to_I_fp_snps, non_A_to_I_fn_snps = set(), set(), set()
gz_flag = False
if hap_vcf.endswith('.gz'):
    gz_flag = True
    f = gzip.open(hap_vcf, 'rb')
else:
    f = open(hap_vcf, 'r')

with f:
    A_to_I_t_tp, A_to_I_t_fp, A_to_I_t_fn = 0, 0, 0
    tp_var_dict, fp_var_dict, fn_var_dict = {}, {}, {}
    tp_var_cnt, fp_var_cnt, fn_var_cnt = 0, 0, 0
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
                tp_var_dict[(ref, alt)] = tp_var_dict.get((ref, alt), 0) + 1
                tp_var_cnt += 1
            else:
                if "FN" in truth and "SNP" in truth:
                    fn_var_dict[(ref, alt)] = fn_var_dict.get((ref, alt), 0) + 1
                    fn_var_cnt += 1
                if "FP" in query and "SNP" in query:
                    fp_var_dict[(ref, alt)] = fp_var_dict.get((ref, alt), 0) + 1
                    fp_var_cnt += 1

            if (ref == 'A' and alt == 'G') or (ref == 'T' and alt == 'C'):
                if "TP" in truth and "SNP" in truth and "TP" in query:
                    A_to_I_t_tp += 1
                else:
                    ## An error genotype will both increase a FP and a FN
                    if "FN" in truth and "SNP" in truth:
                        A_to_I_t_fn += 1
                    if "FP" in query and "SNP" in query:
                        A_to_I_t_fp += 1
            else:
                ## non RNA-editing sites
                if "TP" in truth and "SNP" in truth and "TP" in query:
                    non_A_to_I_tp_snps.add("{}:{}".format(line[0], line[1]))
                    non_A_to_I_tp += 1
                else:
                    ## An error genotype will both increase a FP and a FN
                    if "FN" in truth and "SNP" in truth:
                        non_A_to_I_fn_snps.add("{}:{}".format(line[0], line[1]))
                        non_A_to_I_fn += 1
                        fn_out.write(line_str)
                    if "FP" in query and "SNP" in query:
                        non_A_to_I_fp_snps.add("{}:{}".format(line[0], line[1]))
                        non_A_to_I_fp += 1
                        fp_out.write(line_str)
tp_var_dict = sorted(tp_var_dict.items(), key=lambda x: x[1], reverse=True)
fp_var_dict = sorted(fp_var_dict.items(), key=lambda x: x[1], reverse=True)
fn_var_dict = sorted(fn_var_dict.items(), key=lambda x: x[1], reverse=True)
for var, cnt in tp_var_dict:
    print("TP: {}->{} {:.3f}".format(var[0], var[1], cnt / tp_var_cnt))
print("-" * 20)
for var, cnt in fp_var_dict:
    print("FP: {}->{} {:.3f}".format(var[0], var[1], cnt / fp_var_cnt))
print("-" * 20)
for var, cnt in fn_var_dict:
    print("FN: {}->{} {:.3f}".format(var[0], var[1], cnt / fn_var_cnt))
print("-" * 20)
print(f"A-to-I TP: {A_to_I_t_tp}, A-to-I FN: {A_to_I_t_fn}, A-to-I FP: {A_to_I_t_fp}")
print(f"Non A-to-I TP = {non_A_to_I_tp},\nNon A-to-I FN = {non_A_to_I_fn},\nNon A-to-I FP = {non_A_to_I_fp}")
P = non_A_to_I_tp / (non_A_to_I_tp + non_A_to_I_fp)
R = non_A_to_I_tp / (non_A_to_I_tp + non_A_to_I_fn)
F1 = 2 * P * R / (P + R)
print(f"Recall = {R},\nPrecision = {P},\nF1 = {F1}")
fp_out.close()
fn_out.close()

## evaluate with genotype
gt_errors = (non_A_to_I_fn_snps & non_A_to_I_fp_snps)

tp = len(non_A_to_I_tp_snps)
fp = len(non_A_to_I_fp_snps)
fn = len(non_A_to_I_fn_snps)
print("regular eval:")
print("TP: {}, FN: {}, FP: {}, Recall: {}, Precision: {}, F1: {}".format(tp, fn, fp,
                                                                         tp / (tp + fn),
                                                                         tp / (tp + fp),
                                                                         2 * tp / (2 * tp + fp + fn)))

tp_gt = tp + len(gt_errors)
fp_gt = fp - len(gt_errors)
fn_gt = fn - len(gt_errors)
print("remove genotype eval:")
print("TP: {}, FN: {}, FP: {}, Recall: {}, Precision: {}, F1: {}".format(tp_gt, fn_gt, fp_gt,
                                                                         tp_gt / (tp_gt + fn_gt),
                                                                         tp_gt / (tp_gt + fp_gt),
                                                                         2 * tp_gt / (2 * tp_gt + fp_gt + fn_gt)))

print("genotype accuracy: {}".format(1 - len(gt_errors) / (tp_gt)))


