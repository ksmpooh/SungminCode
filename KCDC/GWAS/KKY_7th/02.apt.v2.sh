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