# 1 Preparing input VCF

## mapping
bwa index ref_genome.fasta -p genome 
ref_genome=/path_to_your_ref_genome/genome
fq_path=/path_to_your_fq/
mkdir bam && cd bam
for file in *.clean.1.fq
do
bname=`basename $file | cut -d. -f1`
bwa mem -t 20 -M $ref_genome -R "@RG\tID:$bname\tPL:ILLUMINA\tSM:$bname" $fq_path/$bname.clean.1.fq $fq_path/$bname.clean.2.fq > $bname.sam && samtools view -bhF 4 -q 30 $i.sam > $bname.F4.bam && samtools sort $bname.F4.bam > $bname.sorted.bam && java -jar picard.jar MarkDuplicates I=$bname.sorted.bam O=$bname.dd.bam M=$bname.dd.metrics REMOVE_DUPLICATES=true  && rm $bname.sam  $bname.F4.bam $bname.sorted.bam &
done

## calling
mkdir vcf && cd vcf
bam_path=/path_to_your_bam/
echo "Chr01 42132932
Chr02	48726069
Chr03	49814079
Chr04	39642072
Chr05	46382547
Chr06	36113639
Chr07	35147422
Chr08	42437421
Chr09	56635340
Chr10	138102" >  chr_info_file.txt
# calling variants from different regions simultaneously
python partition_bcftools_calling.py -c chr_info_file.txt -l 5000000 -b "${bam_path}/*.dd.bam" -t 10 -r "${ref_genome}/genome" > batch_calling.sh
sh batch_calling.sh

## merge
### sort 
vcf_path=/path_to_your_vcf/
for file in $vcf_path*raw.vcf
do
bname=`echo $i | cut -d"_" -f1-3`
echo $bname
bcftools sort ${bname}_raw.vcf -O v -o ${bname}.sort.vcf &
done

### merge by chr
for file in $vcf_path*_r1.sort.vcf
do
i=${file/bcftools_/}
i=${file/_r1.sort.vcf/}
bcftools concat -O v --threads 20 bcftools_${i}_r*.sort.vcf > bcftools_${i}.vcf && bcftools sort bcftools_${i}.vcf -O z -o ${i}.sort.vcf.gz &
done

### merge all
echo "bcftools concat --threads 20 \\" >  mergeVCF.sh
for i in `ls $pwd*.sort.vcf.gz`
do
echo "${i} \\" >> mergeVCF.sh
done 
echo "-O z -o merge.vcf.gz " >>  mergeVCF.sh
sh mergeVCF.sh &

## filtering
### only bi-allelic sites 
### miss rate less than 80
vcftools --gzvcf  merge.vcf.gz  --max-missing 0.80 --minDP 3 --maf 0.05 --mac 3 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out Filtered.millet &