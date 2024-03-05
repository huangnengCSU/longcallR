import argparse


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
                    tp_snps.add("{}:{}".format(line[0], line[1]))
                else:
                    ## An error genotype will both increase a FP and a FN
                    if "FN" in truth and "SNP" in truth:
                        fn_cnt += 1
                        fn_snps.add("{}:{}".format(line[0], line[1]))
                    if "FP" in query and "SNP" in query:
                        fp_cnt += 1
                        fp_snps.add("{}:{}".format(line[0], line[1]))
    return tp_snps, fp_snps, fn_snps


def eval_with_gt(hap_vcf):
    tp_snps, fp_snps, fn_snps = load_hap(hap_vcf)
    gt_errors = (fp_snps & fn_snps)

    tp = len(tp_snps)
    fp = len(fp_snps)
    fn = len(fn_snps)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="hap.py output vcf file", required=True)
    args = parser.parse_args()
    eval_with_gt(args.i)
