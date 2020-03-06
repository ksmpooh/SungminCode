#!/bin/bash

	for chr in $(seq 1 22);do
		plink --bfile /ADATA/smkim/JG/03.QC_2nd/MERGE/JG.merge_rmfrq_rmking --chr $chr --reference-allele /ADATA/smkim/JG/03.QC_2nd/INPUTs/Axiom_KOR.annot.extract.addINDEL.Final.REF.txt --make-bed --out /ADATA/smkim/JG/04.phasing/OUTPUTs/01.split/JG.merge.chr$chr
        done


