import os,glob

wDir = "/DATA/smkim/KKY/05.imputation/"
outDir = wDir + "OUTPUTs/"
shDir = wDir + "SCRIPTs/"
shDir = shDir + "05.vcf.hardy.filter/"
os.system("mkdir %s"%shDir)
outDir = outDir + "5.vcf.merge_MAF0.01INFO0.08/"
hardyDir = outDir + "Hardy_filter/"

os.system("mkdir %s"%hadryDir)


def main():
    dfs = glob.glob(outDir + "*gz")
    for df in dfs:
        vcf = df.replace(outDir,"")
        shout = df.replace(outDir,shDir).replace(".vcf.gz",".hardy.filter.sh")
        out = open(shout,"w")
        hardy_out = df.replace(outDir,hardyDir).replace(".vcf.gz","_hardy")
        out.write("plink --vcf %s --hardy --out %s\n"%(df,hardy_out))
        os.system("awk \'$9<1e-6{print $0}\' test.hwe")

#awk '$9<1e-6{print $0}' test.hwe

#vcftools --gzvcf KKY.7th.20211125.chr22.vcf.gz --hwe 1e-6 --recode --stdout | bgzip -c > test.vcf.gz

#genome@genome109:/DATA/smkim/KKY/05.imputation/OUTPUTs/05.vcf.merge_MAF0.01INFO0.08$ zcat KKY.7th.20211125.chr22.vcf.gz | grep -v "#" |wc -l
111685
zcat test.vcf.gz |grep -v "#" |wc -l
111684


#ls *gz |cut -d"." -f1-4 | xargs -I{} -P 25 bash -c "vcftools --gzvcf {}.vcf.gz --hwe 1e-6 --recode --recode-INFO-all --stdout | bgzip -c > hardy.filter/{}_filter.vcf.gz"

