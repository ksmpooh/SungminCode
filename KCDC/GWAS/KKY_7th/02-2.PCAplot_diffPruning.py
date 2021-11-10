## pruning 방법에 따른 PCA 분석


# python3 /DATA/smkim/KKY/02.1stQC/OUTPUTs/PCA_compare.py input
import os,sys

#ori :   --geno 0.01 --hwe 0.001 --indep-pairwise 50 5 0.02 --maf 0.1 --chr 1-22 --exclude /BDATA/smkim/JG/02.QC_1st/INPUTs/chr6_14_rm.txt
#p1 : --geno 0.1 --hwe 0.001 --maf 0.1 --indep-pairwise 50 5 0.01 --chr 1-22 --exclude /BDATA/smkim/JG/02.QC_1st/INPUTs/chr6_14_rm.txt
#p2 : --geno 0.1 --hwe 0.001 --maf 0.2 --indep-pairwise 50 5 0.01 --chr 1-22 --exclude /BDATA/smkim/JG/02.QC_1st/INPUTs/chr6_14_rm.txt
#p3 : --geno 0.1 --hwe 0.001 --maf 0.4 --indep-pairwise 50 5 0.01 --chr 1-22 --exclude /BDATA/smkim/JG/02.QC_1st/INPUTs/chr6_14_rm.txt

input = sys.argv[1]



def main():
    outDir = "./PCA/"
    out1 = outDir + input + "_t1_pruning"
    out2 = outDir + input + "_t2_pruning"
    out3 = outDir + input + "_t3_pruning"
    os.system("plink --bfile %s --geno 0.1 --hwe 0.001 --maf 0.1 --indep-pairwise 50 5 0.01 --chr 1-22 --exclude /BDATA/smkim/JG/02.QC_1st/INPUTs/chr6_14_rm.txt --out %s"%(input,out1))
    os.system("plink --bfile %s --geno 0.1 --hwe 0.001 --maf 0.2 --indep-pairwise 50 5 0.01 --chr 1-22 --exclude /BDATA/smkim/JG/02.QC_1st/INPUTs/chr6_14_rm.txt --out %s"%(input,out2))
    os.system("plink --bfile %s --geno 0.1 --hwe 0.001 --maf 0.4 --indep-pairwise 50 5 0.01 --chr 1-22 --exclude /BDATA/smkim/JG/02.QC_1st/INPUTs/chr6_14_rm.txt --out %s"%(input,out3))

    os.system("plink --bfile %s --extract %s.prune.in --make-bed --out %s"%(input,out1,out1.replace("pruning","pruned")))
    os.system("plink --bfile %s --extract %s.prune.in --make-bed --out %s"%(input,out2,out2.replace("pruning","pruned")))
    os.system("plink --bfile %s --extract %s.prune.in --make-bed --out %s"%(input,out3,out3.replace("pruning","pruned")))
    #flashpca_x86-64 --bfile KKY.7th.tera_snpolisher_pruned --outpc PCA.txt
    os.system("flashpca_x86-64 --bfile %s --outpc %s"%(out1.replace("pruning","pruned"),outDir+"t1_PCA.txt"))
    os.system("flashpca_x86-64 --bfile %s --outpc %s"%(out2.replace("pruning","pruned"),outDir+"t2_PCA.txt"))
    os.system("flashpca_x86-64 --bfile %s --outpc %s"%(out2.replace("pruning","pruned"),outDir+"t3_PCA.txt"))

main()