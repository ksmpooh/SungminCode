import os,sys,glob


datain = sys.argv[1]


front = 28477797
#6:28477797-33448354     1725    6:1725:A:G
def main():
    vcfin = open(datain,"r")
    vcfout = open(datain.replace(".vcf","_Change.chr.pos.ID.vcf"),"w")

    while 1:
        line = vcfin.readline()
        if line[0] =="#":
            vcfout.write(line)
        else:
            break

    while 1:
        if not line: 
            break
        tmp = line.split("\t")
        Pos = tmp[1]
        tmp[0] = "6"
        tmp[1] = str(front + int(tmp[1]))
        tmp[2] = tmp[2].replace(Pos,tmp[1])
        vcfout.write('\t'.join(tmp))
        line = vcfin.readline()
    
    vcfout.close()

main()

        
    