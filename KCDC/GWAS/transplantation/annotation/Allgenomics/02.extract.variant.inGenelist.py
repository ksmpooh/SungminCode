import os,glob

ref = "/BDATA/smkim/KR_allogenomics/03.transmembrane_uniprot/transmembrane_targetGene.inUniprot.txt"

#1:69610:C:T	missense_variant	OR4F5
#1:738539:T:C	missense_variant	AL669831.1
#1:739132:A:C	missense_variant	AL669831.1
#1:819959:C:T	splice_acceptor_variant	AL645608.2
#1:863258:A:G	missense_variant	AL645608.1


outDir = "gene_extract/"

def main():
    refIn = open(ref,"r")
    refIn = [s.replace("\n","") for s in refIn]
    out= open("tmp.txt","w")
    count = 0
    for line in refIn:
        count = count + 1
        var,t1,gene = line.split()
        if count != 1 and gene != tmp:
            out.close()
            os.system("plink --bfile JG.KR.KD.merge_updateID_extractFuncVariant --extract tmp.txt --recodeA --out ./gene_extract/KR.KD.%s"%(tmp))
            out= open("tmp.txt","w")
        tmp = gene
        out.write("%s\n"%var)


main()
