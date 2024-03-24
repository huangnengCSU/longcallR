import gzip
import sys
def load_callset(input_file, output_file):
    gz_flag = False
    if input_file.endswith(".gz"):
        gz_flag = True
        f = gzip.open(input_file, 'rb')
    else:
        f = open(input_file, 'r')
    fout = open(output_file, 'w')
    for line in f:
        if gz_flag:
            line = line.decode('utf-8')
        if line.startswith("#"):
            fout.write(line)
            continue
        fields = line.strip().split("\t")
        filter = fields[6]
        INFO = fields[7].strip().replace("RDS=", "").split(",")
        if filter == "PASS" and "ase_snp" in INFO:
            fout.write(line)

    f.close()
    fout.close()


if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    load_callset(input_file, output_file)