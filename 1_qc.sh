# 1 Quality control for capture/shotgun metagenomic sequencing for soil samples and WGS sequencing for plant tissues
raw_fq=/path_to_your_raw_data/
clean_fq=/path_to_clean_data/
mkdir -p $clean_fq

# # 1.1 Remove adapters and collapse pair-end reads
cd $raw_fq
for file in $raw_fq/*fastq
do
bname=`basename $file | cut -d. -f1`
AdapterRemoval --file1 $bname.R1.fastq --file2 $bname.R2.fastq --basename $bname --threads 6 --adapter1 AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --trimns --trimqualities --mm 3 --collapse --minalignmentlength 11 --minlength 30 --gzip && zcat $bname.collapsed.gz $bname.collapsed.truncated.gz $bname.pair1.truncated.gz $bname.pair2.truncated.gz $bname.singleton.truncated.gz > $bname.adRemoved.fq && rm $bname.collapsed.gz $bname.collapsed.truncated.gz $bname.discarded.gz $bname.pair1.truncated.gz $bname.pair2.truncated.gz $bname.singleton.truncated.gz && seqkit rename $bname.adRemoved.fq > $bname.adRemoved.1.fq && rm $bname.adRemoved.fq && mv $bname.adRemoved.1.fq $bname.adRemoved.fq && gzip $bname.adRemoved.fq &
done &> adapter_removal_log.txt

# # 1.2 Remove residues and deduplicate
# # # platform: illumina novaseq
echo ">adapter1
AATGATACGGCGACCACCGAGATCTACAC
>adapter2
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>adapter3
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>adapter4
ATCTCGTATGCCGTCTTCTGCTTG
>adapter5
CAAGCAGAAGACGGCATACGAGAT
>adapter6
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>adapter7
AGATCGGAAGAGCGTCGTGTAGGGAAAG AGTGT
>adapter8
GTGTAGATCTCGGTGGTCGCCGTATCATT" > adpater_list.fa

for file in $raw_fq/*fastq
do
bname=`basename $file | cut -d. -f1`
fastp --in1 $bname.adRemoved.fq.gz --out1 $bname.residue_cleaned.fq -w 20 --adapter_fasta adpater_list.fa --html $bname.fastp.html --json $bname.fastp.json --dedup --dup_calc_accuracy 3 --low_complexity_filter --complexity_threshold 30 -l 30 --trim_poly_g --poly_g_min_len 6 --trim_poly_x --poly_x_min_len 6 -p --qualified_quality_phred 30 --average_qual 25 && sga preprocess --dust-threshold=5 -m 30 $bname.residue_cleaned.fq -o $bname.prep.fq && rm $bname.residue_cleaned.fq && vsearch --threads 20 --strand both --minseqlength 30 --fastx_uniques $bname.prep.fq --fastqout $bname.cleaned.fq --log $bname.vsearch_log.txt && rm $bname.prep.fq && gzip $bname.cleaned.fq &
done &> cleaning_log.txt

mv *.cleaned.fq.gz $clean_fq