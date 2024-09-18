# -*- coding:utf-8 -*-
import math
import optparse

usage = """python Partition_Bcftools_calling.py -c chr_info_file -l region_length -b bam_dir -t threads -r reference
                                                                                --Joe"""

parser = optparse.OptionParser(usage)
# 使用 add_option 来定义命令行参数
parser.add_option("-c", dest="chr_info", help="chromosome information file",
                  metavar="FILE", action="store", type="string")
parser.add_option("-l", dest="region_length", help="region length",
                  metavar="STRING", action="store", type="string")
parser.add_option("-b", dest="bam_dir", help="bam directory",
                  metavar="STRING", action="store", type="string")
parser.add_option("-t", dest="threads", help="threads",
                  metavar="STRING", action="store", type="string")
parser.add_option("-r", dest="reference", help="reference",
                  metavar="STRING", action="store", type="string")
(options, args) = parser.parse_args()
# chr_info
# chr_name length
chr_info_file = options.chr_info
region_length = int(options.region_length)
bam_dir = options.bam_dir
threads = options.threads
ref = options.reference

chr_info = {}
with open(chr_info_file) as f:
    for line in f.readlines():
        chr_info[line.strip().split()[0]] = int(line.strip().split()[1])

for chr, length in chr_info.items():
    for i in range(1, math.ceil(length / region_length) + 1):
        if i * region_length < length:
            chr_region = str((i - 1) * region_length + 1) + " - " + str(i * region_length)
        elif i * region_length > length:
            chr_region = str((i - 1) * region_length + 1) + " - " + str(length)
        region = chr + ":" + chr_region
        print(
            "bcftools mpileup -r %s %s  --threads %s --fasta-ref %s | bcftools call -mv  --threads %s -o bcftools_%s_r%s_raw.vcf &" % (
                region, bam_dir, threads, ref, threads, chr, str(i)))