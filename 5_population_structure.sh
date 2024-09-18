# 5 Population structure with captured data

# # 5.1 SNP calling with capture data
ref=/path_to_you_ref/
bam_path=/path_to_your_bam/
vcf_path=/path_to_output_vcf/
mkdir -p $vcf_path
cd $bam_path

# # # index bam files
for file in `ls *dd.bam`
do
samtools index $file
done

mkdir vcf && cd vcf
echo "chr1 length_of_chr1
chr2 length_of_chr2
.
.
chrn length_of_chrn" >  chr_info_file.txt

# # # calling SNP separately
cd $vcf_path
python partition_bcftools_calling.py -c chr_info_file.txt -l 5000000 -b $bam_path/*.dd.bam -t 4 -r $ref > batch_calling.sh
sh batch_calling.sh

# # # merge vcfs 
for file in `ls $pwd*_raw.vcf`
do
bname=`echo $file | cut -d"_" -f1-3`
bcftools sort ${bname}_raw.vcf -O v -o ${bname}.sort.vcf 
done

# # # # merge by chr
for file in `ls $pwd*_r1.sort.vcf`
do
chr=${file/bcftools_/}
chr=${chr/_r1.sort.vcf/}
bcftools concat -O v --threads 20 bcftools_${chr}_r*.sort.vcf > bcftools_${chr}.vcf && bcftools sort bcftools_${chr}.vcf -O z -o ${chr}.sort.vcf.gz &
done

# # # # merge all chrs
echo "bcftools concat --threads 20 \\" >  mergeVCF.sh
for file in `ls $pwd*.sort.vcf.gz`
do
echo "${file} \\" >> mergeVCF.sh
done 
echo "-O z -o captured_merge.vcf.gz " >>  mergeVCF.sh
sh mergeVCF.sh 

# # 5.2 SNP calling with resequencing data
# # # 5.2.1 The mapping and SNP calling of resequencing data is similar to that of captured data except for the collapse in 1_qc.sh

# # 5.3 Extract SNPs in target regions and merge capture and resequencing data
probe_bed=/path_to_your_probe_BED/
bedtools intersect -a captured_merge.vcf.gz -b $probe_bed -sorted -header -wa | bcftools view -Oz -o Cap.vcf.gz
bedtools intersect -a resequencing_merge.vcf.gz -b $probe_bed -sorted -header -wa | bcftools view -Oz -o Res.vcf.gz
tabix -p vcf Cap.vcf.gz
tabix -p vcf Res.vcf.gz
vcf-merge Cap.vcf.gz Res.vcf.gz | bgzip -c > merge_target.vcf.gz

# # 5.4 Filter vcf
vcftools --gzvcf merge_target.vcf.gz --max-missing 0.80 --maf 0.05 --minDP 3 --mac 3 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --plink --out Millet

# # 5.5 IBS distance
plink --vcf $vcf_path/Millet*vcf --make-bed --out Millet
plink --bfile Millet --genome --out ibs_output

# # 5.6 PCA
pop_path=/path_for_population_genetics/
mkdir $pop_path
cd $pop_path
plink --vcf $vcf_path/Millet*vcf --recode --out Millet --const-fid --allow-extra-chr && plink --allow-extra-chr --file Millet --noweb --make-bed --out Millet && plink --allow-extra-chr --threads 36 -bfile Millet --pca 10 --out Millet 

# # 5.7 ADMIXTURE
vcftools --vcf  $vcf_path/Millet*vcf --plink --out Millet
plink --noweb --file Millet --geno 0.1 --maf 0.05 --hwe 0.0001 --make-bed --out Millet --allow-extra-chr
mv Millet.bim old.bim
awk '{$1 = "1"; print $0}' old.bim > Millet.bim
for K in 2 3 4 5 6 7 8
do
nohup admixture --cv Millet.bed $K | tee log${K}.out
done

