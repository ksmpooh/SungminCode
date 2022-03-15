import os,glob

#bcftools view -i 'MAF >= 0.01' JG.KR.NODAT.imputation.chr1.142535439_147000000.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' | awk '{print $1}
##fileformat=VCFv4.1
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
#9	101258881	9_101258881_T_C	T	C	.	.	.
#4	164934874	4_164934874_G_A	G	A	.	.	.
#5	163542505	5_163542505_T_G	T	G	.	.	.

wDir = "/DATA/smkim/JG/anno/"
vcfDir = "/LaCie2/KOTRY/04.imputation/KR.2nd/01.vcf/"
shDir = wDir + "SCRIPTs/"
outDir = wDir + "OUTPUTs/"

#bcftools view -i 'MAF >= 0.01' JG.KR.imputation.mergeGen.processing.chr18.78000001_78017156.vcf.gz |  bcftools query -f '%CHROM\t%POS\t%ID%REF\t%ALT\t.\t.\t.\n'
#JG.KR.imputation.mergeGen.processing.chr18.78000001_78017156.vcf.gz

def main():
    dfs = glob.glob(vcfDir + "*gz")
    for vcf in dfs:
        #vcf = vcf.replace(".tbi","")
        out = vcf.replace(vcfDir,outDir).replace(".vcf.gz","forANNO.vcf")
        os.system("cp %sheader.txt %s"%(wDir,out))
        shout = open(vcf.replace(vcfDir,shDir).replace("vcf.gz","sh"),"w")
        #shout.write("\'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t.\\t.\\t.\\n\'")
        #shout.write("%CHORM")
        shout.write("bcftools view -i \'MAF >= 0.01\' %s | bcftools query -f \'%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t.\\t.\\t.\\n\' >> %s"%(vcf,out))
        #shout.write("bcftools view -i \'MAF >= 0.01\' %s | bcftools query -f \'%%CHROM\\t%%POS\\t%%ID\\"%(vcf))
        #shout.write("bcftools view -i \'MAF >= 0.01\' %s "%(vcf))
        shout.close()

main()
