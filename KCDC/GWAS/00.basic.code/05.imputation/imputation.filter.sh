
#bcftools view -i 'INFO >= 0.8 & MAF >= 0.01'
#tabix -f -p
###bcftools filter
mkdir ./05.vcf.filter_MAF0.01INFO0.8/
mkdir ./05.vcf.filter_INFO0.8/
ls *gz | cut -d"." -f1-5 | xargs -I{} -P 44 bash -c "bcftools view -i 'INFO >= 0.8 & MAF >= 0.01' {}.vcf.gz -Oz > ./05.vcf.filter_MAF0.01INFO0.8/{}_MAF0.01_INFO0.8.filter.vcf.gz"
ls *gz | cut -d"." -f1-4 | xargs -I{} -P 24 bash -c "bcftools view -i 'INFO >= 0.8' {}.vcf.gz -Oz > ./05.vcf.filter_INFO0.8/{}_INFO0.8.filter.vcf.gz"


ls JG*gz | cut -d"." -f1-5 | xargs -I{} -P 24 bash -c "bcftools view -i 'INFO >= 0.8 & MAF >= 0.01' {}.vcf.gz | bcftools sort -Oz > ./05.vcf.filter_MAF0.01INFO0.8/{}_MAF0.01_INFO0.8.filter.vcf.gz"
ls *gz | cut -d"." -f1-4 | xargs -I{} -P 24 bash -c "bcftools view -i 'INFO >= 0.8' {}.vcf.gz | bcftools sort -Oz > ./05.vcf.filter_INFO0.8/{}_INFO0.8.filter.vcf.gz"


ls *gz | cut -d"." -f1-5 | xargs -I{} -P 15  bash -c "bcftools view -i 'R2 >= 0.8' {}.dose.vcf.gz -Oz > ./05.vcf.filter_INFO0.8/{}_INFO0.8.filter.does.vcf.gz"
ls *gz | cut -d"." -f1-4 | xargs -I{} -P 24 bash -c "bcftools view -i 'R2 >= 0.8' {}.vcf.gz -Oz > ./05.vcf.filter_INFO0.8/{}_INFO0.8.filter.vcf.gz"

ls *gz | xargs -I{} -P 4 bash -c 'tabix -f -p vcf {}'
ls *gz | xargs -I{} -P 60 bash -c 'tabix -f -p vcf {}'


 #bcftools sort input.vcf > output.vcf




## ID_ID -> ID 변경 
#bcftools reheader -s [inputfile] -o [new_vcf.gz] [old_vcf.gz]


ls *.gz | cut -d"." -f4 | xargs -I {} -P 22 bash -c "bcftools reheader -h ../new_header_forKKY.6h.txt -o ./FINAL/KBA.KNHANES.6th.Imputed_MINIMAC4.{}.filter_INFO0.8.vcf.gz KKY.6th.imputation_MINIMAC4.{}.filter.vcf.gz"

KBA.KNHANES.6th.Imputed_Minimac4.{}.filter_INFO0.8.vcf.gz
KBA.KNHANES.6th.QCed.PLINK

==================
##VCF sort : bcftools sort input > output
## 먼저 input vcf에서 header를 grep 하여 out.vcf에 저장합니다 
grep "^#" in.vcf > out.vcf 
## 그 다음 header를 제외한 vcf 부분에서 필드 구분자 -k에서 첫 번째 필드 
## -k1에서 1번째 필드를 Version sort 하고, 
## -k2에서 -n numeric sort 해줍니다. 
grep -v "^#" in.vcf | sort -k1,1V -k2n >> out.vcf 
## V 옵션을 넣지 않으면, ## chr10 다음 chr2가 오게 됩니다. 
## n 옵션을 넣지 않으면, ## 숫자가 아닌 문자열로 정렬되어 
## 122134 다음에 12345가 오게 됩니다



#### marker count

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%MAF\t%R2"

cat *txt | awk '$5 >= 0.1 && $6 >= 0.8{print $0}' |wc -l 
cat *txt | awk '$5 >= 0.05 && $5 < 0.1 && $6 >= 0.8{print $0}'|wc -l
cat *txt | awk '$5 >= 0.01 && $5 < 0.05 && $6 >= 0.8{print $0}'|wc -l
cat *txt | awk '$5 < 0.01 && $6 >= 0.8{print $0}'|wc -l


cat *txt | awk '$5 >= 0.1{print $0}' |wc -l 
cat *txt | awk '$5 >= 0.05 && $5 < 0.1{print $0}'|wc -l
cat *txt | awk '$5 >= 0.01 && $5 < 0.05{print $0}'|wc -l
cat *txt | awk '$5 < 0.01{print $0}'|wc -l