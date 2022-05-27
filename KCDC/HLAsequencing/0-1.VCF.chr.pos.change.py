import gzip, os, sys

# python test.py [input.vcf.gz]
# 6:28477797-33448354_104_T_C

"""    2945 6:28477797-33448354     290204  6:28477797-33448354_290204_T_A  T       A       14      .       AF=0.225;AQ=14  GT:DP:AD:GQ:PL:
   2946 6:28477797-33448354     290210  6:28477797-33448354_290210_T_C  T       C       23      .       AF=0.35;AQ=23   GT:DP:AD:GQ:PL:
   2947 6:28477797-33448354     290225  6:28477797-33448354_290225_A_C  A       C       14      .       AF=0.183333;AQ=14       GT:DP:A
   2948 6:28477797-33448354     290237  6:28477797-33448354_290237_C_G  C       G       11      .       AF=0.175;AQ=11  GT:DP:AD:GQ:PL:
   2949 6:28477797-33448354     290255  6:28477797-33448354_290255_G_A  G       A       15      .       AF=0.008333;AQ=15       GT:DP:A
   2950 6:28477797-33448354     290285  6:28477797-33448354_290285_G_A  G       A       29      .       AF=0.333333;AQ=29       GT:DP:A
   2951 6:28477797-33448354     290299  6:28477797 """

vcfIn = sys.argv[1]
print('VCF in : %s'%vcfIn)
vcfOut = vcfIn.replace(".vcf.gz",".updateID.vcf")
def main():
    with open(vcfOut,"w") as out:
        with gzip.open(vcfIn,"rt") as fin:
            for line in fin:
                if line[0] == "#":
                    out.write(line)
                    continue
                pos = line.split("\t")[1]
                new_pos = str(int(pos)+28477797)
                new_line = line.replace("6:28477797-33448354","chr6").replace("\t%s\t"%(pos),"\t%s\t"%(new_pos)).replace("_%s_"%(pos),"_%s_"%(new_pos))
                out.write(new_line)
    
    os.system("bgzip -f %s"%(vcfOut))
    os.system("tabix -f -p vcf %s.gz"%vcfOut)


main()

    



def test():
    wdir = "/Users/ksmpooh/Desktop/KCDC/long_read/2022/VCF/"
    import gzip,os
    count = 0
    with open(wdir + 'test.vcf','w') as fout:
        with gzip.open(wdir + 'HLA.Longread.Seq.Q20.Deepvariant_Variantcalling.GLnexus_Jointcalling.vcf.gz','rt') as fin:
            for line in fin:
                if line[0] == "#":
                    fout.write(line)
                    continue
                count = count + 1
                if count == 10:
                    break
                #print(line[0:200])
                pos = line.split("\t")[1]
                new_line = line.replace("\t%s\t"%(pos),"\t%s\t"%(str(int(pos)+28477797 - 1)))
                fout.write(new_line)
                #print(line.replace("\t104\t","\ttest\t")[0:200])
                #6:28477797-33448354
                
    #os.system("bgzip -c %s)



