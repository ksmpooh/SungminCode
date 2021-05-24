### phasing imputation

import os,sys
wDir = "/DATA/smkim/JG/06.HLAimputation.split/"
inDir = wDir + "INPUTs/"
outDir = wDir + "OUTPUTs/"

def main():
    print("main...:")
    genes = ["A","B","DRB1"]
    genotype_panel = outDir + "01.split/"
    phasing_Dir = outDir + "02.phasing/"
    impute_Dir = outDir + "03.impute4/"
    plink_Dir = outDir + "04.plink/"
    phasing_map = inDir + "map/genetic_map_chr6_combined_b37_addCHR.txt"
    impute4_map = inDir + "map/genetic_map_chr6_combined_b37.txt"
    hap = inDir + "Han.hap.gz"
    legend = inDir + "Han.legend.gz" 
    a1_allele = inDir + "a1.allele.P.ref"
    #JG.QCed.HLA_intersect.bim        JG.QCed.HLA_intersect_HLA.B.bim          JG.QCed.HLA_intersect_HLA.DRB1.bim
    #JG.QCed.HLA_intersect_HLA.A.bim  JG.QCed.HLA_intersect_HLA.B_pruning.bim  JG.QCed.HLA_intersect_HLA.DRB1_pruning.bim

    for gene in genes:
        #phasing_in = genotype_panel + "JG.QCed.HLA_intersect_HLA.%s"%gene
        phasing_in = genotype_panel + "JG.QCed.HLA_intersect_HLA.%s_pruning"%gene 
        phasing_out = phasing_in.replace(genotype_panel,phasing_Dir) 
        #phasing_out = phasing_Dir + "JG.QCed.HLA_intersect_HLA.%s_phasing"%gene
        impute_in = phasing_out + ".haps.gz"
        #impute_out = impute_Dir + "JG.QCed.HLA_intersect_HLA.%s_imputation"%gene
        impute_out = phasing_out.replace(phasing_Dir,impute_Dir)
        #plink_out = plink_Dir + "JG.QCed.HLA_intersect_HLA.%s_imputation"%gene
        plink_out = impute_out.replace(impute_Dir,plink_Dir)
        grep_allele = plink_out + "_onlyTargetAllele.txt"

        os.system("~/Downloads/Eagle_v2.4.1/eagle --bfile %s --geneticMapFile %s --chrom 6 --outPrefix %s --maxMissingPerSnp 0.3 --maxMissingPerIndiv 0.5 --numThreads 8"%(phasing_in,phasing_map,phasing_out))
        os.system("/DATA/smkim/JG/TOOLs/impute4.1.2_r300.2 -no_maf_align -buffer 5000 -int 28477833 33448188 -h %s -l %s -m %s -g %s -o %s"%(hap,legend,impute4_map,impute_in,impute_out))
        os.system("plink --gen %s.gen --sample %s.sample --make-bed --out %s --allow-extra-chr"%(impute_out,phasing_out,plink_out))
        os.system("grep HLA_%s %s.bim | awk '{print $2}' > %s"%(gene,plink_out,grep_allele))
        os.system("plink --bfile %s --a1-allele %s --recodeA --out %s"%(plink_out,a1_allele,plink_out + "_raw"))
#main()




def make_sh():

    theme = "pruning"
    theme = "general"

    print("main...:")
    gene_infos = [["A","28910309","30913647"],["B","30321652","32324956"],["C","31546552","33557625"]]
    #genes = ["A","B","DRB1"]
    genotype_panel = outDir + "01.split/"
    phasing_Dir = outDir + "02.phasing/"
    impute_Dir = outDir + "03.impute4/"
    plink_Dir = outDir + "04.plink/"
    os.system("mkdir " + genotype_panel)
    os.system("mkdir " + phasing_Dir)
    os.system("mkdir " + impute_Dir)
    os.system("mkdir " + plink_Dir)
    phasing_map = inDir + "map/genetic_map_chr6_combined_b37_addCHR.txt"
    impute4_map = inDir + "map/genetic_map_chr6_combined_b37.txt"
    hap = inDir + "Han.hap.gz"
    legend = inDir + "Han.legend.gz" 
    a1_allele = inDir + "a1.allele.P.ref"
    #JG.QCed.HLA_intersect.bim        JG.QCed.HLA_intersect_HLA.B.bim          JG.QCed.HLA_intersect_HLA.DRB1.bim
    #JG.QCed.HLA_intersect_HLA.A.bim  JG.QCed.HLA_intersect_HLA.B_pruning.bim  JG.QCed.HLA_intersect_HLA.DRB1_pruning.bim
    for gene_info in gene_infos:
        gene = gene_info[0]
        front = gene_info[1]
        tail= gene_info[2]
        
        if theme == "pruning":
            phasing_in = genotype_panel + "JG.QCed.HLA_intersect_HLA.%s_pruning"%gene 
        else:
            phasing_in = genotype_panel + "JG.QCed.HLA_intersect_HLA.%s"%gene        
            
        #phasing_in = genotype_panel + "JG.QCed.HLA_intersect_HLA.%s"%gene
        #phasing_in = genotype_panel + "JG.QCed.HLA_intersect_HLA.%s_pruning"%gene 
        phasing_out = phasing_in.replace(genotype_panel,phasing_Dir) 
        #phasing_out = phasing_Dir + "JG.QCed.HLA_intersect_HLA.%s_phasing"%gene
        impute_in = phasing_out + ".haps.gz"
        #impute_out = impute_Dir + "JG.QCed.HLA_intersect_HLA.%s_imputation"%gene
        impute_out = phasing_out.replace(phasing_Dir,impute_Dir)
        #plink_out = plink_Dir + "JG.QCed.HLA_intersect_HLA.%s_imputation"%gene
        plink_out = impute_out.replace(impute_Dir,plink_Dir)
        grep_allele = plink_out + "_onlyTargetAllele.txt"
        
        with open(phasing_in.replace(genotype_panel, "")+".sh","w") as shout:
            shout.write("~/Downloads/Eagle_v2.4.1/eagle --bfile %s --geneticMapFile %s --chrom 6 --outPrefix %s --maxMissingPerSnp 0.3 --maxMissingPerIndiv 0.5 --numThreads 8\n"%(phasing_in,phasing_map,phasing_out))
            shout.write("/DATA/smkim/JG/TOOLs/impute4.1.2_r300.2 -no_maf_align -buffer 5000 -int %s %s -h %s -l %s -m %s -g %s -o %s\n"%(front,tail,hap,legend,impute4_map,impute_in,impute_out))
            shout.write("plink --gen %s.gen --sample %s.sample --make-bed --out %s --allow-extra-chr \n"%(impute_out,phasing_out,plink_out))
            shout.write("grep HLA_%s %s.bim | awk '{print $2}' > %s \n"%(gene,plink_out,grep_allele))
            shout.write("plink --bfile %s --a1-allele %s --extract %s --recodeA --out %s \n"%(plink_out,a1_allele,grep_allele,plink_out + "_raw"))


main_sh()
