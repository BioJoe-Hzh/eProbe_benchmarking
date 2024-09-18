# 2 Taxonomic identification and extract sequences assigned to Setaria genus or within this genus

# # construct the bowtie2 database for mapping
database="/path_to_Poaceae66_database"
clean_fq="path_to_your_clean_data"
out_path="path_to_your_output_dir"
mkdir -p $out_path
names_dmp="/path_to_NCBI_taxonomy/names.dmp"
nodes_dmp="/path_to_NCBI_taxonomy/nodes.dmp"
acc2tax="/path_to_accession2taxid/acc2tax.txt"
thread=

# # make sure that the bowtie2, samtools, and ngsLCA are in your $PATH
for file in $fq_path/*.cleaned.fq.gz
do
sample=`basename $file | cut -d. -f1`
python run_ngsLCA.py -f $sample.cleaned.fq.gz -d $database --thread $thread --names $names_dmp --nodes $nodes_dmp -a $acc2tax --minedit 0 --maxedit 2 -o $out_path/$sample &
done

# # extract sequences assigned to Setaria genus (taxid: 4544) or within this genus
for file in $out_path/*lca
do
sample=`basename $file | cut -d. -f1`
python get_seq_from_lca.py -l $sample.lca -t 4544 -i $fq_path/$sample.cleaned.fq.gz -f
fastq -o $sample.Setaria.fastq &
done