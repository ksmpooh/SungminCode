
#bcftools view -i 'INFO >= 0.8 & MAF >= 0.01'
#tabix -f -p
###bcftools filter
mkdir ./05.vcf.filter_MAF0.01INFO0.08/
mkdir ./05.vcf.filter_INFO0.08/
ls *gz | cut -d"." -f1-4 | xargs -I{} -P 24 bash -c "bcftools view -i 'INFO >= 0.8 & MAF >= 0.01' {}.vcf.gz -Oz > ./05.vcf.filter_MAF0.01INFO0.08/{}_MAF0.01_INFO0.8.filter.vcf.gz"
ls *gz | cut -d"." -f1-4 | xargs -I{} -P 24 bash -c "bcftools view -i 'INFO >= 0.8' {}.vcf.gz -Oz > ./05.vcf.filter_INFO0.08/{}_INFO0.8.filter.vcf.gz"

ls *gz | xargs -I{} -P 24 bash -c 'tabix -f -p vcf {}'


