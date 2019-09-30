# Gastric protocol

## Protocol
<pre><code>

</code></pre>

### 1. Quality control

#### 1.1 1st QC
##### 1.1.1 Genotype calling 
<pre><code>
apt-genotype-axiom --analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis --arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --dual-channel-normalization true --cel-files cel_file_list.txt --summaries --write-models --out-dir outDir/
</code></pre>




##### 1.1.2 SNPolisher
<pre><code>
ps-metrics --posterior-file outDir/AxiomGT1.snp-posteriors.txt --call-file outDir/AxiomGT1.calls.txt --metrics-file outDir/AxiomGT1.out.txt

ps-classification --species-type human --metrics-file outDir/AxiomGT1.out.txt --output-dir outDir/

ps-classification-supplemental --performance-file outDir/Ps.performance.txt --summary-file outDir/AxiomGT1.summary.txt --call-file outDir/AxiomGT1.calls.txt --posterior-file outDir/AxiomGT1.snp-posteriors.txt --output-dir outDir/


</code></pre>

##### 1.1.3  
#### 1.2 2nd QC


<pre><code>

</code></pre>



### 2. Association
#### 2.1 Imputation
#### 2.2 Association
#### 2.3 Annotation & Visualization
## Data handling
## command
## Tool

## Protocol


### 1. Quality control
This process call QC. 
#### 1.1 Genotype calling
To make plink file.
Original chip data format is CEL file. 
We need to change Plink file format(ped and map).
First genotype call using APT ( Affymetrix Power Tool) 

In this step, it needs CEL files and CEL files list txt(colnames : cel_files)


#### 1.2 SNPolisher (by batch)

#### 1.3 Sample QC (by batch)
#### 1.4 merge
#### 1.4 SNP QC

### 2. Association
#### 2.1 Imputation
#### 2.2 Association
#### 2.3 Annotation & Visualization


## Data handling
#### VCF file to data frame in R
VCF (수정필요) vcf file로 실습 안해봄

<pre><code>
a <- read.table("")
#colname
b <- readLines("c:/Users/user/Desktop/snpEff_summary_genome20.genes.txt",n = 2)[2]
c <- gsub("#","",b)   # 파일에 있는 # --> "" 빈공간으로 바꾸기
d <- strsplit(c,"\t") # 다시 그 string을 \t로 구분
e =a
colnames(e)<-unlist(d) #list type으로 되어 있는 d 를 unlist 해준다!
</code></pre>





