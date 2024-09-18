import subprocess
import pandas as pd
import os
import shutil
import datetime


# function: Maps sequences from a FASTA file to a reference database using Bowtie2.
def run_mapping(thread, databse, fastq, output, k):
    """
    Maps sequences from a FASTA file to a reference database using Bowtie2.

    Parameters:
    - thread (int): Number of threads for Bowtie2 to utilize.
    - database (str): Path to the Bowtie2 index (reference database).
    - fasta (str): Path to the input FASTA file containing sequences to map.
    - output (str): Base path for the output SAM file.
    - k (int): The '-k' parameter for Bowtie2, determining the maximum number of alignments to report.
    """
    bname = os.path.basename(databse)
    bowtie2_cmd = f"bowtie2 -k {k} --threads {thread} -x {databse} -U {fastq} --no-unal -S {output}.{bname}.sam"
    subprocess.run(bowtie2_cmd, shell=True, check=True)


# function: Filters and prepares the header for the SAM file to be used with ngsLCA.
def filter_header(header_full, header_new, header_subset):
    """
    Filters and prepares the header for the SAM file to be used with ngsLCA.

    Parameters:
    - header_full (str): Path to the file containing the full header.
    - header_new (str): Path to the file containing new headers to include.
    - header_subset (str): Output path for the filtered header file.
    """
    df1 = pd.read_csv(header_new, header=None)
    df1 = df1[0].unique()
    df1 = [f"SN:{value}" for value in df1]

    header = pd.read_csv(header_full, header=None, sep='\t')
    header = header[header[1].isin(df1)]

    header.to_csv(header_subset, header=False, index=False, sep='\t')


# function: Processes SAM files to prepare for ngsLCA.
def process_sam_files(bname, bDB, thread):
    """
    Processes SAM files to prepare for ngsLCA, including extracting headers, aligning, and converting to BAM.

    Parameters:
    - bname (str): Base name for the output files.
    - bDB (str): Database name used for creating the SAM file.
    - thread (int): Number of threads to use for processing.
    """
    sam_file = f"{bname}.{bDB}.sam"
    header_subset_1_file = f"{bname}.{bDB}.header_subset.1.txt"
    header_subset_2_file = f"{bname}.{bDB}.header_subset.2.txt"
    header_subset_3_file = f"{bname}.{bDB}.header_subset.3.txt"
    alignment_file = f"{bname}.{bDB}.alignment.txt"
    header_new_file = f"{bname}.{bDB}.header_new.txt"
    bam_file = f"{bname}.{bDB}.bam"
    
    # Check if SAM file has any reads mapped (ignoring headers)
    result = subprocess.run(f"samtools view {sam_file} | head -n 1", shell=True, capture_output=True, text=True)
    if not result.stdout.strip():
        print(f"No reads mapped for database {bDB}. Skipping.")
        return None
    
    # Extract header information
    subprocess.run(["grep", "@HD", sam_file], stdout=open(header_subset_1_file, "w"))
    subprocess.run(["grep", "@SQ", sam_file], stdout=open(header_subset_2_file, "w"))
    subprocess.run(["grep", "@PG", sam_file], stdout=open(header_subset_3_file, "w"))
    subprocess.run(["grep", "-v", "^@", sam_file], stdout=open(alignment_file, "w"))

    # Remove duplicate lines from header_subset_2_file
    awk_command = f"awk '!seen[$0]++' {header_subset_2_file} > {header_subset_2_file}.tmp"
    subprocess.run(awk_command, shell=True, check=True)
    os.remove(header_subset_2_file)
    os.rename(header_subset_2_file + ".tmp", header_subset_2_file)

    # Extract unique values from alignment_file
    extract_cmd = f"cat {alignment_file} | cut -f3 | sort -u | uniq > {header_new_file}"
    subprocess.run(extract_cmd, shell=True, check=True)

    # filter header
    os.rename(header_subset_2_file, header_subset_2_file + ".tmp")
    filter_header(header_subset_2_file + ".tmp", header_new_file, header_subset_2_file)

    # remove intermediate files
    # os.remove(header_new_file)
    # os.remove(alignment_file)
    # os.remove(header_subset_2_file + ".tmp")

    # concatenate and convert to BAM
    subprocess.run(
        f"cat {header_subset_1_file} {header_subset_2_file} {header_subset_3_file} {alignment_file} | samtools view -@ {thread} -b -o {bam_file}",
        shell=True, check=True)

    # remove intermediate files
    os.remove(header_subset_1_file)
    os.remove(header_subset_2_file)
    os.remove(header_subset_3_file)
    os.remove(sam_file)


# function: Merges and sorts BAM files, preparing for ngsLCA.
def merge_sort(bam_path, output_base, thread):
    """
    Merges and sorts BAM files, preparing for ngsLCA analysis.

    Parameters:
    - bam_path (str): Directory containing BAM files to merge and sort.
    - bname (str): Base name for the output sorted BAM file.
    - thread (int): Number of threads to use for processing.

    Returns:
    - str: Path to the sorted BAM file.
    """
    # Check the number of BAM files in bam_path
    bam_files = [f for f in os.listdir(bam_path) if f.endswith('.bam')]
    if len(bam_files) == 1:
        # If only one BAM file, skip merging and just sort
        merged_bam = os.path.join(bam_path, bam_files[0])
    else:
        # Merge multiple BAM files
        merge_cmd = f"samtools merge -n -@ {thread} {output_base}.merged {' '.join([f'{bam_path}/{f}' for f in bam_files])}"
        subprocess.run(merge_cmd, shell=True, check=True)
        merged_bam = f"{output_base}.merged"

    # Sort the BAM file
    sort_cmd = f"samtools sort -n -@ {thread} -O bam -o {output_base}.sorted.bam {merged_bam}"
    subprocess.run(sort_cmd, shell=True, check=True)

    return f"{output_base}.sorted.bam"


# function: call ngsLCA to run eDNA pipeline
def run_ngsLCA(bam, names_dmp, nodes_dmp, acc2tax, minedit, maxedit, output):
    """
    Executes ngsLCA for taxonomic identification based on specified parameters.

    Parameters:
    - bam (str): Path to the input BAM file for ngsLCA.
    - names_dmp, nodes_dmp, acc2tax (str): Paths to NCBI taxonomy and accession to taxonomy mapping files.
    - minedit, maxedit (int): Minimum and maximum edit distances for ngsLCA identification.
    - output (str): Base path for the output file generated by ngsLCA.

    Returns:
    - str: Path to the LCA output file.
    """
    ngsLCA_cmd = f"ngsLCA -names {names_dmp} -nodes {nodes_dmp} -acc2tax {acc2tax} -editdistmin {minedit} -editdistmax {maxedit} -bam {bam} -outnames {output}.min{minedit}max{maxedit}"
    subprocess.run(ngsLCA_cmd, shell=True, check=True)
    return f"{output}.min{minedit}max{maxedit}.lca"


# function: Executes another Python script with specified parameters.
def run_script(script_path, parameters):
    """
    Executes another Python script with specified parameters.

    Parameters:
    - script_path (str): Path to the Python script to execute.
    - parameters (list): List of command-line parameters to pass to the script.
    """
    command = ['python', script_path] + parameters
    subprocess.check_output(command, universal_newlines=True)


if __name__ == '__main__':
    import argparse

    usage = """ python Run_ngsLCA.py -f fastq --thread 16 -d ref1,ref2,... --names names.dmp --nodes nodes.dmp -a acc2tax.txt --minedit 0 --maxedit 2 -o output """

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-f", "--fastq", dest="fastq", action="store", nargs='?',
                        help="Input fastq.", metavar="FILE", required=True)
    parser.add_argument("--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Number of thread (default: 1). ", metavar="STRING")
    parser.add_argument("-d", "--database", dest="database", action="store", nargs='?',
                        help="Bowtie2 database(s) (separated with comma).", metavar="STRING", required=True)
    parser.add_argument("-k", dest="keep_hit", action="store", nargs='?', default='100',
                        help="The number of hits to be kept in bowtie2 mapping (default: 100).", metavar="STRING")
    parser.add_argument("--names", dest="names", action="store", nargs='?',
                        help="NCBI taxonomy names.dmp", metavar="FILE", required=True)
    parser.add_argument("--nodes", dest="nodes", action="store", nargs='?',
                        help="NCBI taxonomy nodes.dmp", metavar="FILE", required=True)
    parser.add_argument("-a", "--acc2tax", dest="acc2tax", action="store", nargs='?',
                        help="accession2taxid.txt, the match file of genome accession and taxid.", metavar="FILE",
                        required=True)
    parser.add_argument("--minedit", dest="minedit", action="store", nargs='?', default=0,
                        help="Minimum edit distance to assign a sequence in ngsLCA (default: 0).", metavar="STRING",
                        required=True)
    parser.add_argument("--maxedit", dest="maxedit", action="store", nargs='?', default=2,
                        help="Maximum edit distance to assign a sequence in ngsLCA. All sequences have edit distances with reference genome lower than this value are regarded as the same distance (default: 2).",
                        metavar="STRING", required=True)
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)

    args = parser.parse_args()

    tmp_dir = f"{args.output}_eProbe_temp"
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)

    output_base = os.path.join(tmp_dir, os.path.basename(args.output))

    for db in args.database.split(","):
        print("Mapping probes to database: %s with bowtie2 ..." % db, datetime.datetime.now())
        run_mapping(int(args.thread), db, args.fastq, output_base, args.keep_hit)
    print("Mapping completed ...", datetime.datetime.now())

    for db in args.database.split(","):
        db_name = os.path.basename(db)
        print(f"Preparing input bam file for generated with database: {db} ...", datetime.datetime.now())
        process_sam_files(output_base, db_name, int(args.thread))
    print("Input bam(s) preparation completed ...", datetime.datetime.now())

    print("Merging and sorting bam ...", datetime.datetime.now())
    prepared_bam = merge_sort(tmp_dir, output_base, int(args.thread))
    print("Merging and sorting completed ...", datetime.datetime.now())

    print("Running taxonomic identification with ngsLCA ...", datetime.datetime.now())
    lca_file = run_ngsLCA(prepared_bam, args.names, args.nodes, args.acc2tax, args.minedit, args.maxedit, args.output)
    print("Taxonomic identification completed ...", datetime.datetime.now())
    shutil.rmtree(tmp_dir)
    output_dir = os.path.dirname(args.output)
    files = os.listdir(output_dir)
    for file in files:
        if file.endswith(f"{os.path.basename(args.output)}.sorted.bam.bin"):
            os.remove(file)