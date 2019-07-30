#Genotype call to Plink code

<ped><code>
python DataPrep.py 0 50000 Sample.Info.txt
Axiom_KORV1_1.na35.annot.extaract.txt callDir

</code></ped>

-callDir : Axiom.call.txt directory
-DataPrep.py : 파일 내에 scriptDir 경로 지정

Sample.Info.txt
FAMID	INDID	FAT_ID	MAT_ID	SEX	PHENO
ID001	ID001	0	0	1	1
ID002	ID002	0	0	2	1


