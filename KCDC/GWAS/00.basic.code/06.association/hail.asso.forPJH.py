import hail as hl
hl.init(spark_conf={'spark.driver.memory': '400g'})   # hail 시동 걸고 
# chromosome 범위 설정
chr_list = range(1, 23)
# 결과 파일을 저장할 디렉토리 설정
output_dir = "/RDATA6/LipidGWAS_KBAv1/KoGES/COMMON/HAIL/VCF/LOGI"

# 각 chromosome에 대해 반복
for chr_num in chr_list:
    chr_num = str(chr_num)
    #MT(CHR{chr_num}) = hl.read_matrix_table(f"KBA_CHR{chr_num}.mt")
    MT = hl.read_matrix_table("KBA_CHR"+chr_num+".mt")
    #MT_CHR = hl.read_matrix_table(f"KBA_CHR{chr_num}.mt")
    table = hl.import_table("test.txt", types={"P_ID": hl.tstr, "TUBER": hl.tint32, "SMOKE": hl.tint32, "DRINK": hl.tint32, "HTNMED": hl.tint32, "LIPMED": hl.tint32, "OB": hl.tint32, "SEX": hl.tint32, "AGE": hl.tint32, "COHORT": hl.tint32, "HF": hl.tint32, "NOI": hl.tint32, "DMMED": hl.tint32, "DMINS": hl.tint32, "DMDI": hl.tint32, "DMEX": hl.tint32, "CLIV": hl.tint32, "ORALCON": hl.tint32, "ASTHCU": hl.tint32, "HRTMED": hl.tint32, "MI": hl.tint32, "HT": hl.tint32, "STROKE": hl.tint32, "ISCHEMIC": hl.tint32, "CATA": hl.tint32, "GLAU": hl.tint32, "PDAS": hl.tint32, "PDAL": hl.tint32, "HEMA": hl.tint32, "RENAL": hl.tint32, "ALIV": hl.tint32, "FLIV": hl.tint32, "DM": hl.tint32, "HYLIPID": hl.tint32, "PARKIN": hl.tint32, "ARTH": hl.tint32, "OST": hl.tint32, "GOUT": hl.tint32, "FRAC": hl.tint32, "BIGBABY": hl.tint32, "BRON": hl.tint32, "POL": hl.tint32, "CHRONICSTOM": hl.tint32, "DUOULCER": hl.tint32, "URININFE": hl.tint32, "HYST": hl.tint32, "OVARY": hl.tint32, "GLYCOSURIA": hl.tint32, "TROURIA": hl.tint32, "SMALLBABY": hl.tint32, "PULMON": hl.tint32, "BPH": hl.tint32, "GB": hl.tint32, "PER": hl.tint32, "THYROID": hl.tint32, "GDM": hl.tint32, "GHT": hl.tint32, "ANTICO": hl.tint32, "D1": hl.tint32, "D2": hl.tint32, "PC1": hl.tfloat32, "PC2": hl.tfloat32, "PC3": hl.tfloat32, "PC4": hl.tfloat32, "FID": hl.tstr, "IID": hl.tstr}, missing="NA").key_by("IID")
    
    #MT = MT.annotate_cols(pheno=table[MT_CHR{chr_num}.s])
    MT = MT.annotate_cols(pheno=table[MT.s])

    # logistic regression 실행
    gwas = hl.logistic_regression_rows(y=MT.pheno.RENAL, x=MT.DS, covariates=[1.0, MT.pheno.SEX, MT.pheno.AGE, MT.pheno.PC1, MT.pheno.PC2, MT.pheno.PC3, MT.pheno.PC4], test='wald')
    # 결과를 파일로 내보내기
    #gwas.export(f"{output_dir}/RENAL_chr{chr_num}_out_output.tsv")
    gwas.export(output_dir + "/RENAL_chr"+chr_num+"_out_output.tsv")

