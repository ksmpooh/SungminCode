
# celfiles : /DATA/smkim/KKY/01.GenotypeCalling/INPUTs/cel_files.txt
# output : /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/




##GenotypeCalling
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom 
--analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis 
--arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml 
--dual-channel-normalization true --cel-files /BDATA/smkim/JG/00.rawData/sample_info/KR_2nd_cel_file_list2.txt 
--summaries --write-models --out-dir /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/

##SNPolisher
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics 
--posterior-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.snp-posteriors.txt 
--call-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.calls.txt 
--metrics-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.out.txt


~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification 
--species-type human 
--metrics-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.out.txt 
--output-dir /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental 
--performance-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/Ps.performance.txt 
--summary-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.summary.txt 
--call-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.calls.txt --posterior-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.snp-posteriors.txt 
--output-dir /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/

##v2
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
    libraryDir = "~/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis"
    wDir = "/DATA/smkim/KKY/01.GenotypeCalling/"
    inDir = wDir + "INPUTs/"
    outDir = wDir + "OUTPUTs/"
    celfiles = inDir + "cel_files.txt"
    #callDir = 
    #os.system("%sapt-genotype-axiom --analysis-files-path %s --arg-file %s/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files %s --summaries --write-models --out-dir %s"%(toolDir,libraryDir,libraryDir,celfiles,outDir))
    print("%sapt-genotype-axiom --analysis-files-path %s --arg-file %s/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files %s --summaries --write-models --out-dir %s\n"%(toolDir,libraryDir,libraryDir,celfiles,outDir))
    print("%sps-metrics --posterior-file %sAxiomGT1.snp-posteriors.txt  --call-file %sAxiomGT1.calls.txt  --metrics-file %sAxiomGT1.out.txt\n"%(toolDir,outDir,outDir,outDir))


main()    
