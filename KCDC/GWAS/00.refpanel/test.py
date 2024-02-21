
### split 5K
import os,glob,sys

plink = "KCHIP130K_MERGED_QCed_addINDEL_rmSNP_rmMAF_HLAregion_update"
os.system("mkdir 5Ksplit")
df = open(plink+".fam",'r')
df = df.readlines()

count = 0
g = 0
tmp5K = open("tmp","w")
for i in df:
    count = count + 1
    tmp = i.split()
    tmp5K.write("%s %s\n"%(tmp[0],tmp[1]))
    if count == 5000:
        g = g + 1
        count = 0
        tmp5K.close()
        os.system('plink --bfile %s --keep tmp --make-bed --out 5Ksplit/%s.sample%s'%(plink,plink,str(g)))
        tmp5K = open("tmp","w")


g = g + 1
tmp5K.close()
os.system('plink --bfile %s --keep tmp --make-bed --out 5Ksplit/%s.sample%s'%(plink,plink,str(g)))





