
# celfiles : /DATA/smkim/KKY/01.GenotypeCalling/INPUTs/cel_files.txt
# output : /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/




##GenotypeCalling
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom 
--analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis 
--arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml 
--dual-channel-normalization true --cel-files /DATA/smkim/KKY/01.GenotypeCallingcel_files.txt
--summaries --write-models --out-dir ./

##SNPolisher
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics 
--posterior-file ./AxiomGT1.snp-posteriors.txt 
--call-file ./AxiomGT1.calls.txt 
--metrics-file ./AxiomGT1.out.txt


~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification 
--species-type human 
--metrics-file ./AxiomGT1.out.txt 
--output-dir ./

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental 
--performance-file ./Ps.performance.txt 
--summary-file ./AxiomGT1.summary.txt 
--call-file ./AxiomGT1.calls.txt --posterior-file ./AxiomGT1.snp-posteriors.txt 
--output-dir ./

##v2
/home/genome/Downloads/apt_2.11.3_linux_64_bit_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/ \
--arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
--out-dir ./APT_advnorm/ \
--log-file ./APT_advnorm/apt-genotype-axiom.log \
--allele-summaries true --snp-posteriors-output true --dual-channel-normalization true \
--cel-files ./cel_file_list.txt \
--genotyping-node:snp-priors-input-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.r1.AxiomGT1.models \
--artifact-reduction-trustcheck true \
--artifact-reduction-output-trustcheck true

/home/genome/Downloads/apt_2.11.3_linux_64_bit_x86_binaries/bin/ps-metrics \
--metrics-file ./APT_advnorm/metrics.txt \
--log-file ./APT_advnorm/PS_Metrics.log \
--posterior-file ./APT_advnorm/AxiomGT1.snp-posteriors.txt \
--summary-file ./APT_advnorm/AxiomGT1.summary.txt \
--report-file ./APT_advnorm/AxiomGT1.report.txt \
--call-file ./APT_advnorm/AxiomGT1.calls.txt \
--special-snps /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.r1.specialSNPs

/home/genome/Downloads/apt_2.11.3_linux_64_bit_x86_binaries/bin/ps-classification \
--metrics-file ./APT_advnorm/metrics.txt \
--output-dir ./APT_advnorm/ \
--species-type Human \
--log-file ./APT_advnorm/PS_Class.log \
--ps2snp-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.r1.ps2snp_map.ps \
--priority-order 'PolyHighResolution, NoMinorHom, MonoHighResolution, OTV, CallRateBelowThreshold'


import os,glob

def main():
    toolDir = "~/Downloads/apt_2.11.3_linux_64_bit_x86_binaries/bin/"
    libraryDir = "~/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/"
    wDir = "/DATA/smkim/KKY/01.GenotypeCalling/"
    inDir = wDir + "INPUTs/"
    outDir = wDir + "OUTPUTs/"
    celfiles = inDir + "cel_files.txt"
    #callDir = 
    #os.system("%sapt-genotype-axiom --analysis-files-path %s --arg-file %s/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files %s --summaries --write-models --out-dir %s"%(toolDir,libraryDir,libraryDir,celfiles,outDir))
    print("%sapt-genotype-axiom --analysis-files-path %s \
        --arg-file %sAxiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
            --allele-summaries true --snp-posteriors-output true --dual-channel-normalization true \
                --artifact-reduction-trustcheck true --artifact-reduction-output-trustcheck true \
                --cel-files %s \
                --summaries --write-models \
                    --genotyping-node:snp-priors-input-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.r1.AxiomGT1.models \
                    --out-dir %s\n"%(toolDir,libraryDir,libraryDir,celfiles,outDir))
    #print("%sps-metrics --posterior-file %sAxiomGT1.snp-posteriors.txt  --call-file %sAxiomGT1.calls.txt  --metrics-file %sAxiomGT1.out.txt\n"%(toolDir,outDir,outDir,outDir))
    print("%sps-metrics \
        --metrics-file %smetrics.txt \
        --log-file %sPS_Metrics.log \
        --posterior-file %sAxiomGT1.snp-posteriors.txt \
        --summary-file %sAxiomGT1.summary.txt \
        --report-file %sAxiomGT1.report.txt \
        --call-file %sAxiomGT1.calls.txt \
        --special-snps %sAxiom_KORV1_1.r1.specialSNPs"%(toolDir,outDir,outDir,outDir,outDir,outDir,libraryDir))

main()
