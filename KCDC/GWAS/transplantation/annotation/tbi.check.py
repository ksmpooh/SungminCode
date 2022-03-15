import os,sys,glob

def main():
    tbi = glob.glob("*tbi")
    vcf = glob.glob("*gz")
    out = open("noTBI.list.txt","w")
    for i in vcf:
        if i+".tbi" not in tbi:
            out.write("%s\n"%i)
    out.close()

main()


#bcftools view -i 'MAF >= 0.01' JG.KR.NODAT.imputation.chr1.142535439_147000000.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' | head