TOOLs/qctool  -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.10001_15000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.10001_15000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.15001_20000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.15001_20000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.1_5000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.1_5000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.20001_25000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.20001_25000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.25001_30000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.25001_30000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.30001_35000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.30001_35000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.35001_40000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.35001_40000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.40001_45000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.40001_45000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.45001_50000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.45001_50000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.50001_55000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.50001_55000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.5001_10000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.5001_10000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.55001_60000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.55001_60000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.60001_65000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.60001_65000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.65001_70000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.65001_70000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.70001_75000.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.70001_75000.sampleID.sample -g 05.imputation/OUTPUTs/01.imputation/JG.imputation.chr22.75001_77469.sampleID.51000001_51244237.gen.gz -s 04.phasing/OUTPUTs/03.5Ksplit/JG.phasing.chr22.75001_77469.sampleID.sample -og 05.imputation/OUTPUTs/02.mergeGen/JG.imputation.mergeGen.chr22.51000001_51244237.gen.gz -os 05.imputation/OUTPUTs/02.mergeGen/JG.imputation.mergeGen.chr22.51000001_51244237.sample



zcat 05.imputation/OUTPUTs/02.mergeGen/JG.imputation.mergeGen.chr22.51000001_51244237.gen.gz |cut -d' ' -f 2- |gzip  -c > 05.imputation/OUTPUTs/02.mergeGen/JG.imputation.mergeGen.processing.chr22.51000001_51244237.gen.gz
