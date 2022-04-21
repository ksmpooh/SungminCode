import os, glob


#vcfDir = "/backup/smkim/KKY/6th_20220331/02.imputed.vcf/"
#KNHANES.6th.KBA.V1.1.QCed.fam
#KKY.6th.imputation_MINIMAC4.chr10.filter.vcf.gz
def main():
    vcfs = glob.glob("./*gz")
    #prefix = 
    for vcf in vcfs:
        out = vcf.replace("KKY.6th.imputation_MINIMAC4.","KNHANES.6th.MINIMAC4_Imputed.").replace("filter","filter_INFO0.8")
        os.system("mv %s %s"%(vcf,out))

main()