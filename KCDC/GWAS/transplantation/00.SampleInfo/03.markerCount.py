# SNP count

#bim
#1	AX-100005102	0	69610	T	C


import os,glob

organs = ["KR","KD","LR","LD"]
Data_sets = ["discovery","replication"]
autosomal = [str(a) for a in range(1,23)]
wdir = "/LaCie2/KOTRY/99.open/forOpen/"

out = open("./JG.QCed.SNP.txt","w")
out.write("type\tautosomal\tsex\ttotal\n")
def main():
    for organ in organs:
        for i in Data_sets:
            tmp = "%s%s/%s/01.QCed.PLINK/"%(wdir,organ,i)
            datain = glob.glob(tmp+"*bim")
            print(datain)
            datain = open(datain.pop(),"r")
            auto_count = 0
            sex_count = 0

            while 1:
                line = datain.readline()
                if not line:break
                line = line.split()
                if line[0] in autosomal:
                    auto_count = auto_count + 1
                else:
                    sex_count = sex_count + 1
            
            out.write("%s_%s\t%s\t%s\t%s\n"%(organ,i,str(auto_count),str(sex_count),str(auto_count + sex_count)))


main()
                



            
