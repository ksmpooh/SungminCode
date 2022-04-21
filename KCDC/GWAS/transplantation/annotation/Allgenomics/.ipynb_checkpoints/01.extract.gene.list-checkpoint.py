# extract gene list from uscs ensemble table
# gene chr6 start end
#HLA-DMA chr6 4149562 4169015
#HLA-DMA chr6 4348619 4351570
#HLA-DMA chr6 4149176 4152367
#HLA-DMA chr6 4197855 4199921
#
#


import os, sys

indata = sys.argv[1]
outdata = indata.replace(".txt","_info.txt")

if len(sys.argv) != 3:
	print("python3 python.py info.txt input.bim")
else:
	print("input :"+ sys.argv[1])
	print("ref "+sys.argv[2])


def main():
    df = open(indata,"r")
    out = open(outdata,"w")
    print(out)
    tmp_chrom = ""
    tmp_gene = ""
    tmp_range = []

    while True:
        line = df.readline().replace("\n","")
        if not line:
            break
        line_split = line.split()
        chrom = line_split[1]
        start = line_split[2]
        end = line_split[3]
        gene = line_split[0]
        if tmp_gene != gene:
            if tmp_gene == "":
                out.write("gene chrom start end\n")
                #print(tmp_gene)
            else:
                out.write("%s %s %s %s\n"%(tmp_gene,tmp_chrom,tmp_range[0],tmp_range[1]))
                #print(1)
            tmp_gene = gene
            tmp_chrom = chrom
            tmp_range=[start,end]
        else:
            if start < tmp_range[0]:
                tmp_range[0] = start
            if end > tmp_range[1]:
                tmp_range[1] = end
        #print(line_split[12])
        #print(tmp_range)
        #print(tmp_gene)

    out.close()

main()

indata = outdata
outdata = indata.replace(".txt","_snp.inbim.txt")
refdata = sys.argv[2]
outdata2 = indata.replace(".txt","_snp.inbim_onlySNPID.txt")



def main2():
    df = open(indata,'r')
#    ref = open(refdata,'r')
    out = open(outdata,'w')
    out2 = open(outdata2,"w")
    header = df.readline()
    out.write("gene chrom start end ID pos a1 a2\n")
    while 1:
        line = df.readline().replace("\n","")
        if not line:
            break
        gene,chrom,start,end = line.split()

        ref = open(refdata,'r')
        while 1:
            ref_line = ref.readline().replace("\n","")
            if not ref_line:
                break
            #print(ref_line)
            ref_chrom,ref_ID,ref_zero,ref_pos,ref_a1,ref_a2 = ref_line.split()
            ref_chrom = "chr" + ref_chrom
            if ref_chrom == chrom:
                if int(ref_pos) >= int(start) and int(ref_pos) <= int(end):
                    out.write("%s %s %s %s %s\n"%(line,ref_ID,ref_pos,ref_a1,ref_a2))
                    out2.write("%s\n"%ref_ID)

                else:
                    continue

    out.close()
    out2.close()

main2()