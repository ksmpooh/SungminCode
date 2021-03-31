
cp /BDATA/smkim/imputation.tool.check/INPUTs/KBA/phasing_test/defualt.thread.DS.phasing.test.chr1_10Ksample* /BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/00.haptovcf/

gzip -d /BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/00.haptovcf/defualt.thread.DS.phasing.test.chr1_10Ksample.haps.gz

shapeit -convert \
--input-haps /BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/00.haptovcf/defualt.thread.DS.phasing.test.chr1_10Ksample \
--thread 1 \
--output-vcf /BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/00.haptovcf/defualt.thread.DS.phasing.test.chr1_10Ksample.vcf \
--output-log /BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/00.haptovcf/defualt.thread.DS.phasing.test.chr1_10Ksample.vcf.log