# celfiles : /DATA/smkim/KKY/01.GenotypeCalling/INPUTs/cel_files.txt
# output : /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/




##GenotypeCalling
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom --analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis --arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files /BDATA/smkim/JG/00.rawData/sample_info/KR_2nd_cel_file_list2.txt --summaries --write-models --out-dir /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/

##SNPolisher
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.snp-posteriors.txt --call-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.calls.txt --metrics-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.out.txt

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.out.txt --output-dir /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental --performance-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/Ps.performance.txt --summary-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.summary.txt --call-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.calls.txt --posterior-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.snp-posteriors.txt --output-dir /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/




~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.snp-posteriors.txt --call-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.calls.txt --metrics-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.out.txt

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.out.txt --output-dir /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental --performance-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/Ps.performance.txt --summary-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.summary.txt --call-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.calls.txt --posterior-file /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/AxiomGT1.snp-posteriors.txt --output-dir /BDATA/smkim/JG/01.1stgenocall/OUTPUTs/KR.2nd/



~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file ./AxiomGT1.snp-posteriors.txt --call-file ./AxiomGT1.calls.txt --metrics-file ./AxiomGT1.out.txt 

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file ./AxiomGT1.out.txt --output-dir ./

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental --performance-file ./Ps.performance.txt --summary-file ./AxiomGT1.summary.txt --call-file ./AxiomGT1.calls.txt --posterior-file ./AxiomGT1.snp-posteriors.txt --output-dir ./



## tera
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom --analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis --arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files /DATA/smkim/KKY/01.GenotypeCalling/INPUTs/tera.cel.txt --summaries --write-models --out-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/AxiomGT1.snp-posteriors.txt --call-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/AxiomGT1.calls.txt --metrics-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/AxiomGT1.out.txt 

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/AxiomGT1.out.txt --output-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental --performance-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/Ps.performance.txt --summary-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/AxiomGT1.summary.txt --call-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/AxiomGT1.calls.txt --posterior-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/AxiomGT1.snp-posteriors.txt --output-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyTera/


## DNA
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom --analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis --arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files /DATA/smkim/KKY/01.GenotypeCalling/INPUTs/DNAlink.cel.txt --summaries --write-models --out-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/AxiomGT1.snp-posteriors.txt --call-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/AxiomGT1.calls.txt --metrics-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/AxiomGT1.out.txt 

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/AxiomGT1.out.txt --output-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental --performance-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/Ps.performance.txt --summary-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/AxiomGT1.summary.txt --call-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/AxiomGT1.calls.txt --posterior-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/AxiomGT1.snp-posteriors.txt --output-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19_onlyDNAlink/



## 2nd ALL (rm only low quality samples)
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom --analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis --arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files /DATA/smkim/KKY/01.GenotypeCalling/INPUTs/cel_files_2nd_rmALLmisingHET.txt --summaries --write-models --out-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/AxiomGT1.snp-posteriors.txt --call-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/AxiomGT1.calls.txt --metrics-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/AxiomGT1.out.txt 

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/AxiomGT1.out.txt --output-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental --performance-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/Ps.performance.txt --summary-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/AxiomGT1.summary.txt --call-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/AxiomGT1.calls.txt --posterior-file /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/AxiomGT1.snp-posteriors.txt --output-dir /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.all/



## 2nd DNAlink (rm low quality samples and PCA
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom --analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis --arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files /DATA/smkim/KKY/03.2ndQC/INPUTs/1stQCout/DNAlink/DNAlink.2nd.celfiles.txt --summaries --write-models --out-dir /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/AxiomGT1.snp-posteriors.txt --call-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/AxiomGT1.calls.txt --metrics-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/AxiomGT1.out.txt 

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/AxiomGT1.out.txt --output-dir /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental --performance-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/Ps.performance.txt --summary-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/AxiomGT1.summary.txt --call-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/AxiomGT1.calls.txt --posterior-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/AxiomGT1.snp-posteriors.txt --output-dir /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.DNAlink/

## 2nd Tera (rm low quality samples and PCA)
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom --analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis --arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files /DATA/smkim/KKY/03.2ndQC/INPUTs/1stQCout/Tera/tera.2nd.celfiles.txt --summaries --write-models --out-dir /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/AxiomGT1.calls.txt --metrics-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/AxiomGT1.out.txt 

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/AxiomGT1.out.txt --output-dir /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental --performance-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/Ps.performance.txt --summary-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/AxiomGT1.summary.txt --call-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/AxiomGT1.calls.txt --posterior-file /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/AxiomGT1.snp-posteriors.txt --output-dir /backup/smkim/KKY/01.GenotypeCalling/OUTPUTs/2nd.Tera/



















