#/BDATA/smkim/UK/marker/MAFINFO/IDwithrsID
#/BDATA/smkim/UK/GWAS_cal/marker/onlyID
#/BDATA/smkim/UK/GWAS_cal/marker/trait

import os,glob

#ALT_logz.onlyID.txt
phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]

inDir = "/BDATA/smkim/UK/GWAS_cal/marker/onlyID/"
refDir = "/BDATA/smkim/UK/marker/MAFINFO/IDwithrsID/ukb_mfi_ALL_V3.txt"
shDir = "/BDATA/smkim/UK/SCRIPTs/"
shoutDir = shDir + "merge.withUKrsID/"
os.system("mkdir %s"%shoutDir)
#outDir = "/BDATA/smkim/UK/GWAS_cal/marker/trait/"
rscript = shDir + "marker.merge.withUK.rsID.R"
#Three args : [ref] [input] [output]
def main():
    for pheno in phenoList:
        with open(shoutDir+pheno+"_marker.merge.withUK.rsID.sh","w") as shout:
            shout.write("Rscript --vanilla %s %s %s %s\n"%(rscript, refDir, inDir+pheno+".onlyID.txt",inDir+pheno+"onlyID.withrsID.txt"))

#main()

import os,glob
phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]

inDir = "/BDATA/smkim/UK/130K/marker/onlyID/"
refDir = "/BDATA/smkim/UK/marker/MAFINFO/IDwithrsID/ukb_mfi_ALL_V3.txt"
shDir = "/BDATA/smkim/UK/SCRIPTs/"
shoutDir = shDir + "merge.withUKrsID/"
os.system("mkdir %s"%shoutDir)
#outDir = "/BDATA/smkim/UK/GWAS_cal/marker/trait/"
rscript = shDir + "marker.merge.withUK.rsID.R"

def main2():
    refDir = "/BDATA/smkim/UK/marker/MAFINFO/makeID/"
    outDir = inDir + "splitCHR/"
    os.system("mkdir "+outDir)
    rscript = shDir + "marker.merge.withUK.rsID.V2.R"
    #ukb_mfi_chr16_v3.txt
    #/BDATA/smkim/UK/marker/MAFINFO/makeID/ukb_mfi_chr20_v3.newID.txt
    #UKB_GGT_logz_re.marker.ID_130K.txt
    for pheno in phenoList:
        for chr in range(1,22+1):
            with open(shoutDir+pheno+"_marker.merge.withUK_chr"+str(chr)+".rsID.sh","w") as shout:
                shout.write("Rscript --vanilla %s %s %s %s\n"%(rscript, refDir+"ukb_mfi_chr"+str(chr)+"_v3.newID.txt", inDir+"UKB_"+pheno+"_re.marker.ID_130K.txt",outDir + pheno+".onlyID.withrsID_chr"+str(chr)+".txt"))

main2()