




import os,glob

chunks = ["10178_5000000","5000001_10000000","10000001_15000000","15000001_20000000","20000001_25000000","25000001_30000000","30000001_35000000","35000001_40000000","40000001_45000000","45000001_50000000","50000001_55000000","55000001_60000000","60000001_65000000","65000001_70000000","70000001_75000000","75000001_80000000","80000001_85000000","85000001_90000000","90000001_95000000","95000001_100000000","100000001_105000000","105000001_110000000","110000001_115000000","115000001_120000000","117000001_121485368","142535440_147000000","147000001_152000000","152000001_157000000","157000001_162000000","162000001_167000000","167000001_172000000","172000001_177000000","177000001_182000000","182000001_187000000","187000001_192000000","192000001_197000000","197000001_202000000","202000001_207000000","207000001_212000000","212000001_217000000","217000001_222000000","222000001_227000000","227000001_232000000","232000001_237000000","237000001_242000000","242000001_247000000","245000001_249240543"]

#infoDir = "ID_info/"
#mafDir = "frq/"
#withrefpanelMAF
mafDir = "/BDATA/smkim/imputation.tool.check/INPUTs/gnomad/chunk/"
#chr1_ALL.freq.FINAL.txt
#chr1_75000001_80000000.vcf.gz.frq
wDir = "/BDATA/smkim/imputation.tool.check/OUTPUTs/info_maf/"+panel
panel = "1KGP/"
#outDir = "/BDATA/smkim/imputation.tool.check/OUTPUTs/info_maf/withrefpanelMAF/"
outDir = wDir + "withGnomadMAF/"
os.system("mkdir %s"%outDir)

shDir = wDir + "SCRIPTs/"


#ID	INFO
#1:100000012_G/T	0.99722
#1:100000081_C/T	0.00549
#1:100000185_C/T	0.17395


#CHR                                                                                       SNP   A1   A2          MAF  NCHROBS
#1                                                                            1:94000004:T:C    C    T    0.0003447     5802
#1                                                                        1:94000148:GAATT:G    G GAATT     0.001034     5802

#R CMD BATCH 01.INFO.add.refMAF.R [ref] [input] [output]

def main():
    dfs = glob.glob(wDir + "*txt")
    #ref = mafDir + "chr1_ALL.freq.FINAL.txt"
    ref = mafDir + "gnomAD_MAF_CHR1.ALL.FINAL.txt"
    shoutDir = shDir + "01.Rscript/"
    os.system("mkdir "+shoutDir)
    for df in dfs:
        with open(df.replace(wDir,shoutDir).replace(".txt",".sh"),"w") as shout:
            shout.write("Rscript --vanilla %s01.INFO.add.refMAF.R %s %s %s"%(shDir,ref,df,df.replace(wDir,outDir).replace(".txt",".withMAF.txt")))


main()