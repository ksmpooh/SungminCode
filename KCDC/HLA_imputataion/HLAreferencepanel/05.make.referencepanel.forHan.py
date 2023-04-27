'''
6       31321684        HLA_B*40:03     A       T

6       29910337        SNPS_A_7_29910337_exon1 G       A
6       29910482        SNPS_A_152_29910482_intron1_A   A       T

6       29910385        AA_A_-22_29910338_exon1_V       A       T
6       29910398        AA_A_-2_29910398_exon1_R        A       T
6       29910540        AA_A_3_29910540_exon2   A       T


6       29911991        6:HLA_A_2501:A:P        A       P
6       29911115        6:SNP_A_29911115_G:A:P  A       P
6       29910557        6:SNP_A_29910557:T:A    T       A
6       29910338        6:AA_A_-22_29910338_V:A:P       A       P
6       29910699        6:AA_A_56_29910699:G:R  G       R

'''
#python....

import os,glob

datain = open("test.vcf","r")
dataout = open("test_convert.vcf","w")

while 1:
    line = datain.readline()
    if not line:break
    if "HLA" in line:
        tmp = line.split()
        id = tmp[2].split(":")
        hla,gene,tp = id[1].split("_")
        if len(tp) == 2:
            newid = "HLA_%s*%s"%(gene,tp)
        else:
            td = tp[0:1+1]
            fd = tp[2:]
            newid = "HLA_%s*%s:%s"%(gene,td,fd)
            if newid in ["HLA_DPB1*10:401","HLA_DPB1*10:501","HLA_DPB1*10:601","HLA_DPB1*10:701","HLA_DPB1*13:501","HLA_DPB1*13:801","HLA_DPB1*10:001"]:
                td = tp[0:2+1]
                fd = tp[3:]
                newid = "HLA_%s*%s:%s"%(gene,td,fd)
        dataout.write("%s"%line.replace(tmp[2],newid).replace("\tP\t","\tT\t"))
    elif "SNP" in line:
        #6       29911115        6:SNP_A_29911115_G:A:P  A       P
        #6       29910557        6:SNP_A_29910557:T:A    T       A

        #6       29910337        SNPS_A_7_29910337_exon1 G       A
        #6       29910482        SNPS_A_152_29910482_intron1_A   A       T
        tmp = line.split()
        id = tmp[2].split(":")
        id = id[1].split("_")
        if len(id) == 4:
            newid = "SNPS_%s_%s_%s"%(id[1],id[2],id[3])
            dataout.write("%s"%line.replace(tmp[2],newid).replace("\tP\t","\tT\t"))
        else:
            newid = "SNPS_%s_%s"%(id[1],id[2])
            dataout.write("%s"%line.replace(tmp[2],newid))
    elif "AA_" in line:
#6       29910385        AA_A_-22_29910338_exon1_V       A       T
#6       29910398        AA_A_-2_29910398_exon1_R        A       T
#6       29910540        AA_A_3_29910540_exon2   A       T
#6       29910338        6:AA_A_-22_29910338_V:A:P       A       P
#6       29910699        6:AA_A_56_29910699:G:R  G       R
        tmp = line.split()
        id = tmp[2].split(":")
        id = id[1].split("_")
        if len(id) == 5:
            newid = "AA_%s_%s_%s_%s"%(id[1],id[2],id[3],id[4])
            dataout.write("%s"%line.replace(tmp[2],newid).replace("\tP\t","\tT\t"))
        else:
            newid = "AA_%s_%s_%s"%(id[1],id[2],id[3])
            dataout.write("%s"%line.replace(tmp[2],newid))
    else:
        #6       28480931        6:28480931:C:T  C       T
        #tmp = line.split()
        #chrom,pos,a1,a2 = tmp[2].split(":")
        dataout.write(line.replace(":","_"))

dataout.close()
print("done")




    
