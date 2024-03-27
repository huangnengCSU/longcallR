import sys

import matplotlib.pyplot as plt
import numpy
import pysam


def parse_bam(bam_path):
    errors = []
    qualities = []
    lengths = []
    samfile = pysam.AlignmentFile(bam_path, "rb")
    for read in samfile.fetch():
        if read.is_duplicate or read.is_secondary or read.is_supplementary or read.is_unmapped or not read.has_tag("NM"):
            continue
        query_alignment_length = read.query_alignment_length  ## reference length contains introns, query length is almost the same as lengths of exons
        nm = read.get_tag("NM")
        mapq = read.mapping_quality
        query_length = read.query_length
        error_rate = nm / query_alignment_length * 100
        qualities.append(mapq)
        lengths.append(query_length)
        errors.append(error_rate)
    return errors, qualities, lengths


def plot_error(errors):
    plt.figure()
    plt.hist(errors, range=(0,40), bins=200)
    plt.xlabel("error rate (%)")
    plt.ylabel("number of reads")
    plt.savefig("error_rate.png")


def n50(lengths):
    ## sort contigs longest>shortest
    all_len = sorted(lengths, reverse=True)
    csum = numpy.cumsum(all_len)

    print("N: %d" % int(sum(lengths)))
    n2 = int(sum(lengths) / 2)

    # get index for cumsum >= N/2
    csumn2 = min(csum[csum >= n2])
    ind = numpy.where(csum == csumn2)

    n50 = all_len[ind[0][0]]
    print("N50: %s" % n50)


if __name__ == '__main__':
    bam_path = sys.argv[1]
    errors, qualities, lengths = parse_bam(bam_path)
    plot_error(errors)
    n50(lengths)
