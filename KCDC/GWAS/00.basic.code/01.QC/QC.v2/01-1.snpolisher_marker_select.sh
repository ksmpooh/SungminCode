# snpolisher SNP QC

#1. ori genocall + recommand
#2. ^1 n SSP recommand
#cat <(cat classification/CallRateBelowThreshold.ps classification/Other.ps classification/OffTargetVariant.ps | grep -v probeset_id) <(cat 1st_ssp/classification/Recommended.ps | grep -v probeset_id) | sort |uniq -c | awk '$1==2{print $2}' 
#3. ^(1+2) n adv recommand

#4. plink extract snp


#### bash 01.geno.sh [ori_geno_dir] [ssp_geno_dir] [adv_dir] [merge_output prefix]

if [ $# -ne 4 ]; then
 echo "Usage: $0 [ori_geno_dir] [ssp_geno_dir] [adv_dir]"
 echo "ex) bash 01.geno.sh /DATA/genocall /DATA/genocall/SSP /DATA/genocall/adv KOTRY.AR_2025_KRKD.1stQC" 
 exit -1
else
 echo "ori_geno_dir : $1"
 echo "ssp_geno_dir : $2"
 echo "adv_dir : $3"
 echo "merge_output prefix : $4"
 echo "ok"
fi



ori_geno_dir=$1
ssp_geno_dir=$2
adv_dir=$3


cat <(cat $1/classification/CallRateBelowThreshold.ps $1/classification/Other.ps $1/classification/OffTargetVariant.ps | grep -v probeset_id) <(cat $2/classification/Recommended.ps | grep -v probeset_id) | sort |uniq -c | awk '$1==2{print $2}' > ssp.snp.select.txt
grep -v -f <(cat $1/classification/Recommended.ps $2/classification/Recommended.ps | grep -v probeset_id | sort | uniq -c | awk '{print $2}') <(cat $3/classification/Recommended.ps | grep -v probeset_id) > adv.snp.select.txt

#Axiom_KBAv1.1_SSP
#Axiom_KBAv1.1_adv
#plink --file Axiom_KBAv1.1 --extract ../classification/ --allow-extra-chr --no-fid --no-parents --no-pheno --no-sex --make-bed --out test

plink --file $1/plink/Axiom_KBAv1.1 --extract $1/classification/Recommended.ps --allow-extra-chr --no-fid --no-parents --no-pheno --no-sex --make-bed --out Axiom_KBAv1.1_Original_call
plink --file $2/plink/Axiom_KBAv1.1_SSP --extract ssp.snp.select.txt --allow-extra-chr --no-fid --no-parents --no-pheno --no-sex --make-bed --out Axiom_KBAv1.1_SSP_call
plink --file $3/plink/Axiom_KBAv1.1_adv --extract adv.snp.select.txt --allow-extra-chr --no-fid --no-parents --no-pheno --no-sex --make-bed --out Axiom_KBAv1.1_adv_call

echo "Axiom_KBAv1.1_SSP_call" > merge.list
echo "Axiom_KBAv1.1_adv_call" >> merge.list

plink --bfile Axiom_KBAv1.1_Original_call --merge-list merge.list --make-bed --out $4



