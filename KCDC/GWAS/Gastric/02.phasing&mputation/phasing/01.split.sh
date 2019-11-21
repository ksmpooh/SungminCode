#!/bin/bash

	for chr in $(seq 1 22);do
		plink --bfile /DATA/smkim/Gastric/last_qc_2nd/OUTPUTs/MERGE/gastric.merge_rmking --chr $chr --reference-allele /DATA/smkim/Gastric/last_qc_2nd/INPUTs/Axiom_KOR.annot.extract.addINDEL.Final.REF.txt --make-bed --out /DATA/smkim/Gastric/Phasing/OUTPUTs/splitPlink/gastric.merge.chr$chr
        done

