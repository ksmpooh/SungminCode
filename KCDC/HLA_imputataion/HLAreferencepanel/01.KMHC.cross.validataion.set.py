import os,glob,sys

## KMHC croos validation set
## 5 fold
## HLA-TAPAS

### 1. nomen (data split)
### 2. make reference
### 3. SNP2HLA imputation
outDir = "/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/KMHC_5CV"
shDir = "/BDATA/smkim/JG.HLAimputation/HLA-TAPAS"



#nomen_data = '/BDATA/smkim/JG.HLAimputation/HLA-TAPAS/DATA/01.nomenclean/HLA.4digit.529sample.nomenclean.chped'
nomen_data = '/BDATA/smkim/JG.HLAimputation/HLA-TAPAS/DATA/01.nomenclean/HLA.4digit.520sample.nomenclean.chped'
plink_data = '/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/00.rawData/JG.forHLAtyping.MHCref.28477797_33448354'
theme = plink_data.replace("/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/00.rawData/JG.forHLAtyping.","")

index_data = '/BDATA/smkim/JG.HLAimputation/HLA-TAPAS/DATA/00.rawDATA/Final_520sample.index.txt'

index_len = range(1,5+1)
'''index_data
ID index
ID1 1
ID2 2
ID3 3
'''

## data split with R

def data_pro():
    print("index,data")
    for i in index_len:
        ref_in = open(nomen_data,'r')
        ref = ref_in.readlines()
        print("awk '$2==%s{print $0,$0}' %s > tmp"%(str(i),index_data))
        os.system("awk '$2==%s{print $1,$1}' %s > tmp"%(str(i),index_data))
        tmp = open("tmp","r")
        tmp = tmp.readlines()
        sampleList = [s.replace("\n","").split()[0] for s in tmp]
        out_path = '%s/01.nomen/g%s_HLA.4digit.nomenclean.chped'%(outDir,str(i))
        out = open(out_path,"w")
        for j in ref:
            tmp = j.split()[0]
            if tmp not in sampleList:
                out.write(j)
            else:
                pass

        out.close()
        os.system("plink --bfile %s --remove tmp --make-bed --out %s/00.KBA/g%s_HLAimputation_train.%s"%(plink_data,outDir,str(i),theme))
        os.system("plink --bfile %s --keep tmp --make-bed --out %s/00.KBA/g%s_HLAimputation_test.%s"%(plink_data,outDir,str(i),theme))
        
'''
python -m MakeReference \
    --variants  MakeReference/example/g1k_subset_snps\
    --chped MakeReference/example/g1k_subset.chped \
    --hg 19 \
    --out MakeReference/example/g1k_subset.bglv4 \
    --dict-AA MakeReference/data/hg19/HLA_DICTIONARY_AA.hg19.imgt3320 \
    --dict-SNPS MakeReference/data/hg19/HLA_DICTIONARY_SNPS.hg19.imgt3320 \
    --phasing
'''
def makeReference(i):
    v_in = outDir + "/00.KBA/g%s_HLAimputation_train.%s"%(str(i),theme)
    chped_in = outDir + "/01.nomen/g%s_HLA.4digit.nomenclean.chped"%str(i)
    out = outDir + "/02.makeReference/g%s_HLAreference.Panel.KBA.%s"%(str(i),theme)
    os.system("mkdir %s"%(outDir + "/02.makeReference"))
    return "python -m MakeReference --variants %s --chped %s --hg 19 --out %s --dict-AA MakeReference/data/hg19/HLA_DICTIONARY_AA.hg19.imgt3320 --dict-SNPS MakeReference/data/hg19/HLA_DICTIONARY_SNPS.hg19.imgt3320 --phasing --nthreads 8 --mem 16g\n"%(v_in,chped_in,out)
'''
python -m SNP2HLA \
    --target SNP2HLA/example/1958BC \
    --reference resources/1000G.bglv4 \
    --out MySNP2HLA/IMPUTED.1958BC \
    --nthreads 2 \
    --mem 4g
'''
def SNP2HLA(i):
    v_in = outDir + "/00.KBA/g%s_HLAimputation_test.%s"%(str(i),theme)
    ref_in = outDir + "/02.makeReference/g%s_HLAreference.Panel.KBA.%s"%(str(i),theme)
    out = outDir + "/03.SNP2HLA/g%s_HLAimputation.SNP2HLA.KBA.%s"%(str(i),theme)
    os.system("mkdir %s"%(outDir + "/03.SNP2HLA"))
    return "python -m SNP2HLA --target %s --reference %s --out %s --nthreads 8 --mem 16g\n"%(v_in,ref_in,out)




def main():
    print("main")
    #index_dataIn = open(index_data,"r")
    #index_dataIn = index_dataIn.readlines()
    data_pro()
    
    for i in index_len:
        out = open("%s/g%s.sh"%(shDir,str(i)),'w')
        out.write(makeReference(i))
        out.write(SNP2HLA(i))
        out.close()

    #DATA split

main()




'''
def data_pro():
    print("index,data")
    # index 별 train/test set 만들기
    for i in range(1,5+1):
        # data 불러오기
        ref_in = open("HLA.4digit.520sample.nomenclean.chped",'r')
        ref = ref_in.readlines()
        
        # index 별 샘플 정보 추출
        print("awk '$2==%s{print $0,$0}' %s > tmp"%(str(i),index_data))
        os.system("awk '$2==%s{print $1,$1}' %s > tmp"%(str(i),index_data))

        # 추출한 정보 읽기
        tmp = open("tmp","r")
        tmp = tmp.readlines()
        sampleList = [s.replace("\n","").split()[0] for s in tmp]

        # HLA type data 나누기
        out_path = '%s/01.nomen/g%s_HLA.4digit.nomenclean.chped'%(outDir,str(i))
        out = open(out_path,"w")
        for j in ref:
            tmp = j.split()[0]
            if tmp not in sampleList:
                out.write(j)
            else:
                pass
        out.close()

        # PLINK data 나누기
        os.system("plink --bfile %s --remove tmp --make-bed --out %s/00.KBA/g%s_HLAimputation_train.%s"\
            %(plink_data,outDir,str(i),theme))
        os.system("plink --bfile %s --keep tmp --make-bed --out %s/00.KBA/g%s_HLAimputation_test.%s"\
            %(plink_data,outDir,str(i),theme))
'''        
