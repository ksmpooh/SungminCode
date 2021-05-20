import os

infoDir = "/BDATA/myhwang/UK/MAFINFO/info0.8/" #Info0.8_ukb_mfi_chr17_v3.txt
outDir = "/backup/smkim/UK/"
bgenDir ="/BDATA/myhwang/UK/IMP/" #ukb_imp_chr19_v3.bgen
shDir = "/BDATA/smkim/UK/SCRIPTs/bgen.filter/"
os.system("mkdir %s"%shDir)

tool = "/BDATA/myhwang/UK/IMP/test/qctool"
#./qctool -g ukb_imp_chr22_v3.bgen -og qctool.filterUsingrsID.test.bgen -incl-rsids rsID.list -threads 16

def main():
	for i in range(1,22+1):
		info = infoDir + "Info0.8_ukb_mfi_chr%s_v3.txt"%str(i)
		bgen = bgenDir + "ukb_imp_chr%s_v3.bgen"%str(i)
		out = outDir + "ukb_imp_chr%s_v3_info0.8.bgen"%str(i)
		with open(shDir + "bgen.filter.chr%s.sh"%(str(i)),"w") as sh:
			sh.write("%s -g %s -incl-rsids %s -og %s"%(tool,bgen,info,out))

main()


def bgen_filter(inDir,outDir,refDir,shDir,Tool):
    print("bgen_filter....")
    for i in range(1,22+1):
        shout = open(shDir + "UKB.bgen.filter.chr%s.sh"%str(i))



#./qctool -g ukb_imp_chr22_v3.bgen -og qctool.filterUsingrsID.test.bgen -incl-rsids rsID.list -threads 16
def sh():
    server = "omics"
    #server = "109"
    if server == "omics":
        bgenDir = ""
        outDir = ""
        refDir = ""
        shDir = ""
        Tool = ""
    else:
        bgenDir = ""
        outDir = ""
        refDir = ""
        shDir = ""
        Tool = ""