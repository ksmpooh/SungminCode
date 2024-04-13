import os,glob,sys
## new final 
## KMHC croos validation set
## 5 fold
## HLA-TAPAS

### 1. PLINK data split
### 2. eaple phasing (ref Guide)
### 3. minimac4 imputation

outDir = "/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/KMHC_IMGT3320_notrare/03.HLAimputation/"
phasingDir = outDir + "01.phasing/"
imputationDir = outDir + "02.imputation/"
refDir = outDir + "99.ref/"
preimpDir = outDir + "00.preimp/"

os.system("mkdir %s %s %s %s %s"%(outDir,phasingDir,imputationDir,refDir,preimpDir))
shDir = "/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/KMHC_IMGT3320_notrare/SCRIPTs/"

#v1 = '/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/00.rawData/new/JG.forHLAtyping.v2023_snpolisher_rmaffy_indel_flip_rmdup_idchange_rmSampleQCout_25M-35M_updateids_HLAtypesample_snpQC'
#v2 = '/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/00.rawData/new/JG.forHLAtyping.v2023_snpolisher_rmaffy_indel_flip_rmdup_idchange_rmSampleQCout_25M-35M_updateids_snpQC_HLAtypesample'
version = "v1"
#version = "v2"

if version == 'v1':
    plink_data = '/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/00.rawData/new/JG.forHLAtyping.v2023_snpolisher_rmaffy_indel_flip_rmdup_idchange_rmSampleQCout_25M-35M_updateids_snpQC_HLAtypesample'
    ref_all = '/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/KMHC_IMGT3320_notrare/02.makeReference/notrare/KMHC_v1.hg19.HLAtapas_makereference_2field_notrare.bgl.phased.vcf.gz'
    theme = "SNPQC_520sample"
elif version == 'v2':
    plink_data = '/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/00.rawData/new/JG.forHLAtyping.v2023_snpolisher_rmaffy_indel_flip_rmdup_idchange_rmSampleQCout_25M-35M_updateids_HLAtypesample_snpQC'
    ref_all = '/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/KMHC_IMGT3320_notrare/02.makeReference/notrare/KMHC_v2.hg19.HLAtapas_makereference_2field_notrare.bgl.phased.vcf.gz'
    theme = "520sample_SNPQC"
else:
    print("Error!!!!!")
    exit()

#theme = plink_data.replace("/BDATA/smkim/JG.HLAimputation/HLAreferencePanel/00.rawData/JG.forHLAtyping.","")
index_data = '/BDATA/smkim/JG.HLAimputation/HLA-TAPAS/DATA/00.rawDATA/Final_520sample.index.txt'

index_len = range(1,5+1)
'''index_data
ID index
ID1 1
ID2 2
ID3 3
'''

## data split with R

def data_pro(i):
    print("index,data")

    print("awk '$2==%s{print $0,$0}' %s > tmp"%(str(i),index_data))
    os.system("awk '$2==%s{print $1,$1}' %s > tmp"%(str(i),index_data))
    os.system("awk '$2==%s{print $1}' %s > tmp1"%(str(i),index_data))
    plink_out = "%s/g%s_HLAimputation_test.%s"%(preimpDir,str(i),theme)
    os.system("plink --bfile %s --keep tmp --make-bed --out %s"%(plink_data,plink_out))
    os.system("bgzip -c %s.vcf > %s.vcf.gz"%(plink_out,plink_out))
    os.system("tabix -f -p vcf %s.vcf.gz"%(plink_out))
    os.system("plink --bfile %s --recode vcf --out %s"%(plink_out,plink_out))
    ref_out = "%s/g%s_HLAimputation_train.%s.vcf.gz"%(preimpDir,str(i),theme)
    os.system("bcftools view -S ^tmp1 %s -Oz -o %s"%(ref_all,ref_out))
    os.system("tabix -f -p vcf %s"%ref_out)
    return ref_out,plink_out


'''
os.system(tool + " --vcfRef "+ref+" --vcfTarget "+i+".vcf.gz --geneticMapFile "+m+" --chrom "+chr+" --vcfOutFormat z --outPrefix /SDATA/smkim/KBA_130K/04.phasing/03.refguide.chr_phasing/KCHIP130K_phasing_MERGED_QCed_addINDEL_rmSNP_rmMAF_chr" +chr+" --numThreads 90")
'''
def phasing(i,ref_out,plink_out):
    m = "/BDATA/smkim/GWAS/ref/map/genetic_map_chr6_combined_b37_addCHR.txt"
    phasing_out = plink_out.replace(preimpDir,phasingDir).replace("_test","_test_phasing")
    return "~/Downloads/Eagle_v2.4.1/eagle --vcfRef "+ref_out+" --vcfTarget "+plink_out+".vcf.gz --geneticMapFile "+m+" --chrom 6 --vcfOutFormat z --outPrefix "+phasing_out.replace(".vcf","")+" --numThreads 20\n"

def imputation(i,phasing_out,ref_out):
    imputation_out = phasing_out.replace(phasingDir,imputationDir).replace("_test_phasing",".minimac4.KBA").replace(".vcf.gz","")
    return "/BDATA/smkim/TOOLs/Minimac4 --minRatio 0.000001 --rsid \
            --refhaps "+ref_out+".m3vcf.gz --haps "+phasing_out+".vcf.gz --noPhoneHome --allTypedSites \
            --format GT,DS,GP --prefix "+imputation_out+" --mapFile /BDATA/smkim/GWAS/ref/map/m3vcf/genetic_map_chr6_combined_b37_addCHR.m3vcf.txt --referenceEstimates --cpu 15\n"

'''
/BDATA/smkim/TOOLs/Minimac4 --minRatio 0.000001 --window 500000 \
--refhaps /SDATA/smkim/KBA_130K/12.panel/m3vcfs/chr6_wgs8k_imputationPanel.m3vcf.gz --haps KCHIP130K_phasing_MERGED_QCed_addINDEL_rmSNP_rmMAF_chr6.vcf.gz --noPhoneHome --allTypedSites \
--format GT,DS,GP --prefix ../../05.imputation/KCHIP130K_imputation_MERGED_QCed_addINDEL_rmSNP_rmMAF_chr6.vcf.gz --mapFile /BDATA/smkim/GWAS/ref/map/m3vcf/genetic_map_chr6_combined_b37_addCHR.m3vcf.txt --referenceEstimates --cpu 15
'''


def main():
    print("main")
    #index_dataIn = open(index_data,"r")
    #index_dataIn = index_dataIn.readlines()
    
    for i in index_len:
        ref_out,plink_out = data_pro(i)
        out = open("%s/g%s.sh"%(shDir,str(i)),'w')
        out.write(phasing(i,ref_out,plink_out))
        phasing_out = plink_out.replace(preimpDir,phasingDir).replace("_test","_test_phasing")
        out.write("/BDATA/smkim/TOOLs/Minimac3/bin/Minimac3 --refHaps %s --rsid --processReference --prefix %s\n"%(ref_out,ref_out))
        out.write(imputation(i,phasing_out,ref_out))
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
