import hail as hl
hl.init(spark_conf={'spark.driver.memory': '400g'})


inDir = "/ADATA/smkim/KBA_130K/05.imputation/"
table = (hl.import_table('/ADATA/smkim/KBA_130K/00.input/transformation_lipid_kchip130k_v1v2_20200414_forhail.ped',impute=True).key_by("IND_ID_IND_ID"))

for i in range(6,6+1):
    inVCF = inDir + "KCHIP130K_imputation_MERGED_QCed_addINDEL_rmSNP_rmMAF_chr%s.vcf.gz.dose.vcf.gz"%str(i)
    inhailMT = inDir + 'hail/KMCHP130K_8Kimp_chr%s_INFO0.8.mt'%str(i)
    #hl.import_vcf(inVCF,force_bgz=True).write(inhailMT, overwrite=True)
    mt = hl.import_vcf(inVCF,force_bgz=True)
    mt = mt.filter_rows(mt.row.info.R2 >= 0.8)
    mt.write(inhailMT, overwrite=True)
    
    print("read Hail Matrix:chr %s\n"%(str(i)))
    mt = hl.read_matrix_table(inhailMT)
    
    print("Annotate table:chr %s\n"%(str(i)))
    mt = mt.annotate_cols(pheno=table[mt.s])
    
    print("Ass TC z")
    outAsso = "/ADATA/smkim/KBA_130K/06.asso/KCHIP130K_8Kimp_Asso_TC_z_chr%s.txt"%str(i)
    gwas = hl.linear_regression_rows(y=mt.pheno.TC_z,x=mt.DS,covariates=[1.0])
    gwas.export(outAsso)
    
    print("Ass HDL z")
    outAsso = "/ADATA/smkim/KBA_130K/06.asso/KCHIP130K_8Kimp_Asso_HDL_z_chr%s.txt"%str(i)
    gwas = hl.linear_regression_rows(y=mt.pheno.HDL_z,x=mt.DS,covariates=[1.0])
    gwas.export(outAsso)
    
    print("Ass LDL z")
    outAsso = "/ADATA/smkim/KBA_130K/06.asso/KCHIP130K_8Kimp_Asso_LDL_z_chr%s.txt"%str(i)
    gwas = hl.linear_regression_rows(y=mt.pheno.LDL_z,x=mt.DS,covariates=[1.0])
    gwas.export(outAsso)
    
    print("Ass TC logz")
    outAsso = "/ADATA/smkim/KBA_130K/06.asso/KCHIP130K_8Kimp_Asso_TG_logz_chr%s.txt"%str(i)
    gwas = hl.linear_regression_rows(y=mt.pheno.TG_logz,x=mt.DS,covariates=[1.0])
    gwas.export(outAsso)