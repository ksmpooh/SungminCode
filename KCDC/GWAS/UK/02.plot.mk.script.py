#### plot

#os.system("Rscript --vanilla 03.ploting.R %s %s"%(datain,plotout))
#/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/assoMERGE_new/merge
# UKB_MAP.nRES_all.MERGE.for.manhattan.txt
import os,glob

phenoList = ["GLU_inv", "HbA1c_inv", "DBP.nRES", "SBP.nRES","PP.nRES","MAP.nRES","HT"]


def main1():
    for pheno in phenoList:
        #datain = "UKB_%s_all.MERGE.for.manhattan.txt"%pheno
        datain = "UKB_%s_ALL_MERGE_filINFO0.8_MAF0.01.for.manhattan.txt"%pheno
        plotout = datain.replace(".txt","")
        cmd = "Rscript --vanilla 03.ploting.R %s %s"%(datain,plotout)
        with open(plotout + ".sh","w") as out:
            out.write(cmd)

main1()

def main2():
    for pheno in phenoList:
        datain = "UKB_%s_all.MERGE.for.manhattan.txt"%pheno
        plotout = datain.replace(".txt","")
        os.system("Rscript --vanilla 03.ploting.R %s %s"%(datain,plotout))
        
#main()