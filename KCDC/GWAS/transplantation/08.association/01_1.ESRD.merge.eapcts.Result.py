###

import os,glob
wDir = "/BDATA/smkim/JG/08.asso/OUTPUTs/"

epactsDir = wDir + "01.b.firth/"
mergeDir = wDir + "01_1.merge/"



def main():
    phenos = ["HT","HT_T2D","sub_Total","T2D","Total"]
    #for pheno in phenos:
    for pheno in ["HT_T2D"]:
        print("pheno : " + pheno)
        inDir = epactsDir + pheno + "/"
        output = mergeDir + "ESRD.%s.association.epacts.merge.txt"%pheno
        os.system("cp %sheader.txt %s"%(mergeDir,output))
        os.system("zcat %s*.epacts.gz | grep -v ^# >> %s"%(inDir,output))

main()

