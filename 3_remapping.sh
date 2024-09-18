# 3 re-mapping and filtering
ref=/path_to_your_reference_genome/
fq_path=/path_to_your_Setaria_fq/
bam=/path_to_your_bam/
mkdir -p $bam

cd $bam
for file in $fq_path/*.Setaria.fastq
do
bname=`basename $file | cut -d. -f1`
bowtie2 -x $ref -p 16 --end-to-end -S $bam/$bname.sam -U $fq_path/$bname.Setaria.fastq && samtools view -F 4 -Sb -q 30 -@ 6 $bam/$bname.sam > $bam/$bname.bam && samtools sort $bam/$bname.bam > $bam/$bname.sorted.bam && java -jar ~/software/picard.jar MarkDuplicates I=$bam/$bname.sorted.bam O=$bam/$bname.dd.bam M=$bam/$bname.dd.metrics REMOVE_DUPLICATES=true && rm  $bam/$bname.sam $bam/$bname.bam $bam/$bname.sorted.bam &
done
