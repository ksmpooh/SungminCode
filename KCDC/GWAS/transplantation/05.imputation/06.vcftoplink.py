### vcf to plink

import os,sys,glob

def main():
    vcfDir = "/LaCie2/KOTRY/04.imputation/KR.2nd/01.vcf/"
    outDir = "/ADATA/smkim/JG/05.vcftoplink/OUTPUTs/"
    dfs = glob.glob(vcfDir + "*vcf.gz")
    for df in dfs:
        with open(df.replace(vcfDir,"").replace(".vcf.gz",".sh"),'w') as shout:
            shout.write('plink --vcf %s --make-bed --out %s'%(df,df.replace(vcfDir,outDir).replace(".vcf.gz","")))


#main()
# bed, bim,fam
def plink_merge():
    wDir = "/BDATA/smkim/JG/05.vcftoplink/"
    plinkDir = wDir + "OUTPUTs/"
    outDir = plinkDir + "02.chrMerge/"
    os.system("mkdir "+outDir)
    inDir = wDir + "INPUTs/"
    os.system("mkdir "+inDir)
    for chr in range(3,22+1):
        dfs = glob.glob("/ADATA/smkim/JG/05.vcftoplink/OUTPUTs/*chr%s.*.bim"%chr)
        i = dfs.pop()
        inout_path = inDir + "merge_list_chr%s.txt"%(str(chr))
        with open(inout_path,"w") as inout:
            for df in dfs:
                inout.write("%s\t%s\t%s\n"%(df.replace(".bim",".bed"),df,df.replace(".bim",".fam")))
        out = outDir + "JG.KR.imputation.plink.chr%s"%str(chr)
        os.system("plink --bfile %s --merge-list %s --allow-no-sex --make-bed --out %s"%(i.replace(".bim",""),inout_path,out))

plink_merge()        
        
    


