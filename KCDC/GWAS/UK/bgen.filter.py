import os

#infoDir = "/BDATA/myhwang/UK/MAFINFO/info0.8/" #Info0.8_ukb_mfi_chr17_v3.txt
#outDir = "/backup/smkim/UK/"
#bgenDir ="/BDATA/myhwang/UK/IMP/" #ukb_imp_chr19_v3.bgen
#shDir = "/BDATA/smkim/UK/SCRIPTs/bgen.filter/"
#os.system("mkdir %s"%shDir)

#tool = "/BDATA/myhwang/UK/IMP/test/qctool"
#./qctool -g ukb_imp_chr22_v3.bgen -og qctool.filterUsingrsID.test.bgen -incl-rsids rsID.list -threads 16

def info_filter():
	for i in range(1,22+1):
		info = infoDir + "Info0.8_ukb_mfi_chr%s_v3.txt"%str(i)
		bgen = bgenDir + "ukb_imp_chr%s_v3.bgen"%str(i)
		out = outDir + "ukb_imp_chr%s_v3_info0.8.bgen"%str(i)
		with open(shDir + "bgen.filter.chr%s.sh"%(str(i)),"w") as sh:
			sh.write("%s -g %s -incl-rsids %s -og %s"%(tool,bgen,info,out))


def filein(datain):
    a = open(datain,'r')
    return [s.replace("\n","") for s in a]

def bgen_filter(inDir,outDir,refDir,shDir,Tool):
    print("bgen_filter....")
    for i in range(1,22+1):
        ref = filein(refDir + "UKB_CHR%s_CHUNK.txt"%str(i))
        for j in ref:
            chr,front,tail = j.split("_")
            bgenIn = inDir + "ukb_imp_%s_v3.bgen"%chr
            bgenOut = outDir + "ukb_imp_%s_v3.bgen"%j
            with open(shDir + "UKB.bgen.filter.%s.sh"%j,'w') as shout:
                shout.write("%s -g %s -incl-range %s-%s -og %s"%(Tool,bgenIn,front,tail,bgenOut))
                


#ukb_imp_chr1_v3.bgen
#head UKB_CHR14_CHUNK.txt
#chr14_19000017_20037108
#chr14_20037153_20636854
#chr14_20636870_21161211
#chr14_21161217_21698759
#chr14_21698837_22269835
#chr14_22269868_22796612
#chr14_22796626_23346317
#qctool -g ../ukb_imp_chr22_v3_info0.8.bgen -incl-range 10000000-20000000 -og chr22_split_test.bgen
#./qctool -g ukb_imp_chr22_v3.bgen -og qctool.filterUsingrsID.test.bgen -incl-rsids rsID.list -threads 16

def sh():
    #server = "omics"
    #server = "109"
    server = "102"
    #server = "OAS"
    if server == "omics":
        bgenDir = ""
        outDir = ""
        refDir = ""
        shDir = ""
        Tool = ""
    elif server == "OAS":
        bgenDir = "/jdata/scratch/myhwang/UK/IMP/"
        outDir = "/jdata/scratch/myhwang/UK/IMP_filter/"
        refDir = "/jdata/scratch/myhwang/KBA_130K/11_UKB/INPUTs/"
        shDir = ""
        Tool = "/jdata/scratch/myhwang/TOOLs/qctool_v2.0-rc9-CentOS6.8-x86_64/qctool"
    else:
        bgenDir = "/BDATA/myhwang/UK/"
        outDir = "/BDATA/myhwang/UK/IMP_filter/"
        refDir = "/BDATA/myhwang/KBA_130K/11_UKB/INPUTs/"
        shDir = "/BDATA/smkim/UK/SCRIPTs/bgen.filter/"
        Tool = "/BDATA/smkim/JG/TOOLs/qctool"
	#Tool = "/jdata/scratch/myhwang/TOOLs/qctool_v2.0-rc9-Ubuntu16.04-x86_64/qcttol"
    bgen_filter(bgenDir, outDir, refDir, shDir, Tool)

sh()



