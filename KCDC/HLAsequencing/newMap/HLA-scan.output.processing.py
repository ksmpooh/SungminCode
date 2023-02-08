##### HLA-scan output 정리
#####
# python HLA-scan.output.processing.py [sampleID] [HLA-scan.output.print.file]
# python HLA-scan.processing.py NIH19KT0248 HLA.scan.NIH19KT0248.txt
# ls *txt |  cut -d"." -f3 | xargs -I {} -P 1 bash -c "python3 HLA-scan.processing.py {} HLA.scan.{}.txt"
import os,sys,glob

def main():
    ID = sys.argv[1]
    filein = sys.argv[2]
    fileout = filein.replace(".txt",".result_processing.txt")
    print("Sample ID : %s\nFile In : %s\nFile Out : %s"%(ID,filein,fileout))
    datain = open(filein,'r')
    df = datain.readlines()
    df.append('HLAscan End')
    
    ### 변수초기화
    out = []
    header = "ID"
    gene = ""
    allele = []
    HLAscan_count = 0
    tmp = ""
    for i in df:
        if "HLAscan" in i:
            if HLAscan_count == 0: #시작
                HLAscan_count = HLAscan_count + 1
                continue
            else:
                header = header + "\t%s_1\t%s_2"%(gene,gene)
                if len(allele) == 0:
                    tmp = tmp + "\tNA\tNA"
                elif len(allele) == 1:
                    tmp = tmp + "\t%s\tNA"%allele[0]
                else:
                    tmp = tmp + "\t%s\t%s"%(allele[0],allele[1])
            if "End" in i:
                print(i)
                print(tmp)
                out.append(header + "\n")
                out.append(tmp + "\n")
                break
            allele_count = 0
            gene = ""
            allele = []
            continue
        if "NIH" in i:
            ID = i.replace("\n","")
            print("ID : %s"%ID)
            #tmp = tmp + "\t%s"%ID
            tmp = "%s"%ID
            continue
        if 'HLA gene' in i:
            #print(i.replace("\n",""))
            print(i.split(":")[1].strip())
            gene = i.split(":")[1].strip()
            continue
        if "[Type 1" in i:
            #print(i.replace("\n",""))
            print(i.split()[2])
            #tmp = tmp + "\t" + i.split()[2]
            allele.append(i.split()[2])
            continue
        if "[Type 2" in i:
            #print(i.replace("\n",""))
            print(i.split()[2])
            allele.append(i.split()[2])
            continue
    print("Done")

    with open(fileout,'w') as dataout:
        for i in out:
            dataout.write(i)
    print("Write Done...")

main()






