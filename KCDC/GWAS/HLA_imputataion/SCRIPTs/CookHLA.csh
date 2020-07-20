#!/bin/csh -f

#################################################################################################################################
#
# COOKHLA: Imputation of HLA amino acids and classical alleles from SNP genotypes
#
# Author: Seungho Cook (kukshomr@gmail.com)
#       
## USAGE: "USAGE: ./COOKHLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink  geneticMap {optional: java_max_memory[mb]}"
#
#
#################################################################################################################################
#
#Thank you for downloading COOKHLA. To use this package:
#
# INPUTS: 
# 1. Plink dataset
# 2. Reference dataset (beagle format)
# 
# DEPENDENCIES: (download and place in the same folder as this script)
# 1. PLINK (1.07)
# 2. Beagle (4.1)
# 3. merge_tables.pl (Perl script to merge files indexed by a specific column)
# 4. linkage2beagle and beagle2linkage (Beagle utilities for PED <-> Beagle format)
# 5. beagle2vcf and vcf2beagle (Beagle utilities for Beagle format <-> vcf)
# 6. excluding_snp_and_refine_target_position-v1COOK02222017.R
# 7. redefineBPv1BH.py
# 8. Panel-subset.py
# 9. bgl2GC_trick_bgl-v1.2BH-07052917.py and excluding_snp_and_refine_target_position-v1COOK02222017.R and GC_tricked_bgl2ori_bgl-v1.2BH-07052917.py,complete_header.R
# 10.Panel-BGL2BED.sh

## USAGE: "USAGE: ./COOKHLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink  geneticMap average_erate {optional: java_max_memory[mb]}"
#
######################################################################################################



if ($#argv < 6) then
    echo "USAGE: ./SNP2HLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink  geneticMap average_erate {optional: java_max_memory[mb]}"; exit 1
endif



set SCRIPTPATH=`dirname $0`

set MERGE=$SCRIPTPATH/merge_tables.pl
set PARSEDOSAGE=$SCRIPTPATH/ParseDosage.csh
set refine=$SCRIPTPATH/redefineBPv1BH.py
set Panelsubset=$SCRIPTPATH/Panel-subset.py
set BGL2BED=$SCRIPTPATH/Panel-BGL2BED.sh
set bgl2GC_trick_bgl=$SCRIPTPATH/bgl2GC_trick_bgl-v1.2BH-07052917.py
set GC_trick_bgl2ori_bgl=$SCRIPTPATH/GC_tricked_bgl2ori_bgl-v1.2BH-07052917.py
# CHECK FOR DEPENDENCIES
if (! -e `which $4`) then
    echo "Please install PLINK (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) and point to the plink run file."; 
    echo "tcsh: use plink"
    echo "bash: use ./plink"
    exit 1
else if (! -e $SCRIPTPATH/beagle4.jar) then
    echo "Please install Beagle 4.1 (http://faculty.washington.edu/browning/beagle/beagle.html#download) and copy the run file (beagle.3.0.4/beagle.jar) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/linkage2beagle.jar) then
    echo "Please copy linkage2beagle.jar included in the beagle 3.0.4 zip file (beagle.3.0.4/utility/linkage2beagle.jar) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/beagle2linkage.jar) then # We use beagle2linkage (Buhm, 8/13/12)
    echo "Please copy beagle2linkage.jar (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) into $SCRIPTPATH/"; exit 1
else if (! -e $MERGE) then
    echo "Please copy merge_tables.pl (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $PARSEDOSAGE) then
    echo "Please copy ParseDosage.csh (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $refine) then
    echo "Please copy ParseDosage.csh (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $Panelsubset) then
    echo "Please copy Panelsubset.py (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/beagle2vcf.jar) then
    echo "Please copy beagle2vcf.jar (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/vcf2beagle.jar) then
    echo "Please copy vcf2beagle.jar (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/bgl2GC_trick_bgl-v1.2BH-07052917.py) then
    echo "Please copy bgl2GC_trick_bgl-v1.2BH-07052917.py (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/excluding_snp_and_refine_target_position-v1COOK02222017.R) then
    echo "Please copy excluding_snp_and_refine_target_position-v1COOK02222017.R (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/GC_tricked_bgl2ori_bgl-v1.2BH-07052917.py) then
    echo "Please copy GC_tricked_bgl2ori_bgl-v1.2BH-07052917.py(included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/complete_header.R) then
    echo "Please copy complete_header.R(included with this package) into $SCRIPTPATH/"; exit 1   
else if (! -e $BGL2BED) then
    echo "Please copy Panel-BGL2BED.sh(included with this package) into $SCRIPTPATH/"; exit 1       
     
endif

# INPUTS
set INPUT=$1
set REFERENCE=$2
set OUTPUT=$3
set PLINK=$4
set geneticMap=$5
set AVER_ERATE=$6

if ($#argv >= 7) then
    set MEM=$7
else
    set MEM=2000 # Default java memory 2000 Mb (2Gb)
endif


set JAVATMP=$OUTPUT.javatmpdir
mkdir -p $JAVATMP
alias plink '$PLINK --noweb --silent --allow-no-sex'
alias beagle 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/beagle4.jar'
alias linkage2beagle 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/linkage2beagle.jar'
alias beagle2linkage 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/beagle2linkage.jar'
alias beagle2vcf 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/beagle2vcf.jar'
alias vcf2beagle 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/vcf2beagle.jar'
alias excluding_target_snp_not_reference 'Rscript $SCRIPTPATH/excluding_snp_and_refine_target_position-v1COOK02222017.R'
alias complete_header.R 'Rscript $SCRIPTPATH/complete_header.R'

# Functions to run
set EXTRACT_MHC = 1
set FLIP        = 1
set CONVERT_IN  = 1
set IMPUTE      = 1
set CONVERT_OUT = 1
set CLEANUP     = 1

# SET PARAMETERS
set TOLERATED_DIFF = .15
set i = 1

echo ""
echo "SNP2HLA: Performing HLA imputation for dataset $INPUT";
echo "- Java memory = "$MEM"Mb"


set MHC=$OUTPUT.MHC

if ($EXTRACT_MHC) then
    echo "[$i] Extracting SNPs from the MHC."; @ i++
    plink --bfile $INPUT --chr 6 --from-mb 29 --to-mb 34 --maf 0.025 --make-bed --out $OUTPUT.MHC
endif
	
if ($FLIP) then
    echo "[$i] Performing SNP quality control."; @ i++

    # Identifying non-A/T non-C/G SNPs to flip
    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.bim >> $OUTPUT.tmp1
    echo "SNP 	POSR	A1R	A2R" > $OUTPUT.tmp2
    cut -f2,4- $REFERENCE.bim >> $OUTPUT.tmp2
    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP |  grep -v -w NA > $OUTPUT.SNPS.alleles

    awk '{if ($3 != $6 && $3 != $7){print $1}}' $OUTPUT.SNPS.alleles > $OUTPUT.SNPS.toflip1
    plink --bfile $MHC --flip $OUTPUT.SNPS.toflip1 --make-bed --out $MHC.FLP

    # Calculating allele frequencies
    plink --bfile $MHC.FLP --freq --out $MHC.FLP.FRQ
    sed 's/A1/A1I/g' $MHC.FLP.FRQ.frq | sed 's/A2/A2I/g' | sed 's/MAF/MAF_I/g' > $OUTPUT.tmp

    mv $OUTPUT.tmp $MHC.FLP.FRQ
    $MERGE $REFERENCE.FRQ.frq $MHC.FLP.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.frq
    sed 's/ /\t/g' $OUTPUT.SNPS.frq | awk '{if ($3 != $8){print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" 1-$10 "\t*"}else{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 "\t."}}' > $OUTPUT.SNPS.frq.parsed
    
    # Finding A/T and C/G SNPs
    awk '{if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")){if ($4 > $7){diff=$4 - $7; if ($4 > 1-$7){corrected=$4-(1-$7)}else{corrected=(1-$7)-$4}}else{diff=$7-$4;if($7 > (1-$4)){corrected=$7-(1-$4)}else{corrected=(1-$4)-$7}};print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" diff "\t" corrected}}' $OUTPUT.SNPS.frq.parsed > $OUTPUT.SNPS.ATCG.frq

    # Identifying A/T and C/G SNPs to flip or remove
    awk '{if ($10 < $9 && $10 < .15){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toflip2
    awk '{if ($4 > 0.4){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toremove

    # Identifying non A/T and non C/G SNPs to remove
    awk '{if (!(($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C"))){if ($4 > $7){diff=$4 - $7;}else{diff=$7-$4}; if (diff > '$TOLERATED_DIFF'){print $1}}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    #awk '{if (($2 != "A" && $2 != "C" && $2 != "G" && $2 != "T") || ($3 != "A" && $3 != "C" && $3 != "G" && $3 != "T")){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 == $5 && $3 != $6) || ($3 == $6 && $2 != $5)){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove

    # Making QCd SNP file
    plink --bfile $MHC.FLP --geno 0.2 --exclude $OUTPUT.SNPS.toremove --flip $OUTPUT.SNPS.toflip2 --make-bed --out $MHC.QC
    plink --bfile $MHC.QC --freq --out $MHC.QC.FRQ
    sed 's/A1/A1I/g' $MHC.QC.FRQ.frq | sed 's/A2/A2I/g' | sed 's/MAF/MAF_I/g' > $OUTPUT.tmp
    mv $OUTPUT.tmp $MHC.QC.FRQ.frq
    $MERGE $REFERENCE.FRQ.frq $MHC.QC.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.QC.frq

    cut -f2 $OUTPUT.SNPS.QC.frq | awk '{if (NR > 1){print $1}}' > $OUTPUT.SNPS.toinclude

    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.QC.bim >> $OUTPUT.tmp1

    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP | awk '{if (NR > 1){if ($5 != "NA"){pos=$5}else{pos=$2}; print "6\t" $1 "\t0\t" pos "\t" $3 "\t" $4}}' > $MHC.QC.bim

    # Recoding QC'd file as ped
    plink --bfile $MHC.QC --extract $OUTPUT.SNPS.toinclude --make-bed --out $MHC.QC.reorder
    plink --bfile $MHC.QC.reorder --recode --out $MHC.QC

    # Making SNP files
    awk '{print "M " $2}' $MHC.QC.map > $MHC.QC.dat
    cut -f2 $MHC.QC.map > $MHC.snps
    cut -d ' ' -f1-5,7- $MHC.QC.ped > $MHC.QC.nopheno.ped

    # Remove temporary files
	rm $OUTPUT.tmp1 $OUTPUT.tmp2
	rm $MHC.FLP*
	rm $MHC.QC.ped $MHC.QC.map
	rm $OUTPUT.SNPS.*
endif

if ($CONVERT_IN) then
    echo "[$i] Convering data to beagle format."; @ i++
    linkage2beagle pedigree=$MHC.QC.nopheno.ped data=$MHC.QC.dat beagle=$MHC.QC.bgl standard=true > $OUTPUT.bgl.log

    echo "===============================================================================" >> $OUTPUT.bgl.log

	#rm $MHC.QC.reorder*
	rm $MHC.QC.nopheno.ped
    
    echo "[$i] Convering data to reference_markers_Position"; @ i++
    
    $refine $REFERENCE.markers $REFERENCE.refined.markers
    
    echo "[$i] Convering data to target_markers_Position and extract not_including snp"; @ i++
    
    #awk '{print $2" "$4" "$5" "$6}' $MHC.QC.bim > $MHC.QC.markers
    awk '{print $2" "$4" "$5" "$6}' $MHC.QC.reorder.bim > $MHC.QC.markers
    rm $MHC.QC.reorder*
    
    excluding_target_snp_not_reference $MHC.QC.markers $REFERENCE.refined.markers $MHC.QC.pre.markers
    
    mv $MHC.QC.bgl $MHC.QC.pre.bgl.phased
    
    awk '{print $1}' $MHC.QC.pre.markers > selected_snp.txt
    
    $Panelsubset $MHC.QC.pre all selected_snp.txt $MHC.QC.refined
    
    
    echo "[$i] Convering data to GC_change_beagle format"; @ i++
    
    #target
    
    $bgl2GC_trick_bgl $MHC.QC.refined.bgl.phased $MHC.QC.refined.markers $MHC.QC.GCchange.bgl $MHC.QC.GCchange.markers
    
    #reference
    
    $bgl2GC_trick_bgl $REFERENCE.bgl.phased $REFERENCE.refined.markers $REFERENCE.GCchange.bgl.phased $REFERENCE.GCchange.markers
    
    
    echo "[$i] Convering data to vcf_format"; @ i++
    
    #target
    beagle2vcf 6 $MHC.QC.GCchange.markers $MHC.QC.GCchange.bgl 0 > $MHC.QC.vcf
    
    #reference
    beagle2vcf 6 $REFERENCE.GCchange.markers $REFERENCE.GCchange.bgl.phased 0 > $REFERENCE.vcf
    
    echo "[$i] Convering data to reference_phased"; @ i++
    
    sed "s%/%|%g" $REFERENCE.vcf > $REFERENCE.phased.vcf
    
    awk '{print $1" "$2" "$3}' $geneticMap > $geneticMap.first
    awk '{print $2}' $REFERENCE.GCchange.markers > $geneticMap.second
    paste -d " " $geneticMap.first $geneticMap.second > $geneticMap.refined.map
    
	rm $geneticMap.first
	rm $geneticMap.second
    
    
    
endif

if ($IMPUTE) then
    echo "[$i] Performing HLA imputation (see $OUTPUT.MHC.QC.imputation_out.log for progress)."; @ i++
    
    set aver_erate = `cat $AVER_ERATE`
	
	
    
    beagle gt=$MHC.QC.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.imputation_out impute=true gprobs=true lowmem=true map=$geneticMap.refined.map ne=10000 overlap=5000 err=$aver_erate
    gzip -d $MHC.QC.imputation_out.vcf.gz
    
    
    
    
endif


if ($CONVERT_OUT) then

	echo "[$i] Converting imputation vcf to beagle"; @ i++
	
	cat $MHC.QC.imputation_out.vcf | vcf2beagle 0 $MHC.QC.imputation_GCchange
	
	gzip -d $MHC.QC.imputation_GCchange.bgl.gz
	
	echo "[$i] Converting imputation GC_beagle to ori_beagle"; @ i++
	
	$GC_trick_bgl2ori_bgl $MHC.QC.imputation_GCchange.bgl $REFERENCE.refined.markers $MHC.QC.imputation_ori.bgl
	
	complete_header.R $MHC.QC.GCchange.bgl $MHC.QC.imputation_ori.bgl HLA_IMPUTED_Result.$MHC.bgl.phased
	
	echo "[$i] Converting imputation genotypes to PLINK .ped format."; @ i++
	
    cat HLA_IMPUTED_Result.$MHC.bgl.phased | beagle2linkage $OUTPUT.tmp # Buhm
    cut -d ' ' -f6- $OUTPUT.tmp.ped > $OUTPUT.tmp       # Buhm
    paste -d ' ' $MHC.fam $OUTPUT.tmp | tr -d "\015" > HLA_IMPUTED_Result.$MHC.ped
    cut -f1-4 $REFERENCE.bim > HLA_IMPUTED_Result.$MHC.map
    cp $MHC.fam HLA_IMPUTED_Result.$MHC.fam

    # Create PLINK bed file
    plink --ped HLA_IMPUTED_Result.$MHC.ped --map HLA_IMPUTED_Result.$MHC.map --make-bed --out HLA_IMPUTED_Result.$MHC
	
	rm HLA_IMPUTED_Result.$MHC.bim
	
	cp $REFERENCE.bim HLA_IMPUTED_Result.$MHC.bim

endif


if ($CLEANUP) then

	echo "[$i] CLEANUP";
	rm -r $OUTPUT.javatmpdir
	rm $geneticMap.refined.map
	rm $REFERENCE.vcf
	rm $REFERENCE.refined.markers
	rm $REFERENCE.phased.vcf
	rm $REFERENCE.GCchange.bgl.phased
	rm $REFERENCE.GCchange.markers
	rm $OUTPUT.* 
#	rm -r $JAVATMP
	rm -f plink.log
	rm selected_snp.txt
    echo "DONE!"
    echo ""
endif
