'''
V1                          V2                          V3
1 HLA_A*01:01 NIH19KT0009_NIH19KT0009=0|0 NIH19KT0003_NIH19KT0003=0|0
2 HLA_A*02:01 NIH19KT0009_NIH19KT0009=0|0 NIH19KT0003_NIH19KT0003=1|0
3 HLA_A*02:06 NIH19KT0009_NIH19KT0009=0|0 NIH19KT0003_NIH19KT0003=0|0
4 HLA_A*02:07 NIH19KT0009_NIH19KT0009=0|0 NIH19KT0003_NIH19KT0003=0|0
5 HLA_A*03:01 NIH19KT0009_NIH19KT0009=0|0 NIH19KT0003_NIH19KT0003=0|0
'''

import os,sys

# python 04.check.GT.py input.txt output.txt
fileIn = sys.argv[1]
refIn = "/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/Final/F1_score/HLAtype.KMHC.GT.after.txt"
ref = open(refIn,"r")
ref = ref.readlines()
print(fileIn)
out = fileIn.replace(".txt",".after.txt")
out = open(out,'w')
def data_pro():
    #header = "ID\t"
    df = open(fileIn,"r")
    
    df = df.readlines()
    header = "\t".join([s.strip().split("_")[0] for s in df[0].split()])
    samples = [s.strip().split("_")[0] for s in df[0].split()[1:]]
    hla_types = []
    out.write(header+"\n")
    for i in df:
        tmp = i.split()
        new_line = tmp[0]
        hla_types.append(tmp[0])
        for j in tmp[1:]:
            tmp2 = j.split("=")[1]
            if tmp2 == "0|1":
                new_line= new_line + "\t1|0"
            else:
                new_line= new_line + "\t%s"%tmp2
        out.write("%s\n"%new_line)
    
    return samples,hla_types

samples,hla_types =  data_pro()
out.close()
newIn = fileIn.replace(".txt",".after.txt")
new_out = newIn.replace(".txt",".GTmatching.txt")
out = open(new_out,"w")
infer = open(newIn,"r")
infer = infer.readlines()

def mk_dic(tmp_df):
    new_dic = {}
    header = tmp_df[0]
    header = header.split()
    for i in tmp_df[1:]:
        i = i.strip().split()
        new_dic[i[0]] = {}
        for id,gt in zip(header[1:],i[1:]):
            new_dic[i[0]][id] = gt
    return new_dic

def check_GT(real_GT,infer_GT):
    #print(real_GT,infer_GT)
    TN = 0
    TP = 0
    FN = 0
    FP = 0
    if real_GT == "0|0":
        if infer_GT == "0|0":
            TN = 2
        elif infer_GT in ["1|0","0|1"]:
            TN = 1
            FP = 1
        elif infer_GT == "1|1":
            FP = 2
    elif real_GT in ["1|0","0|1"]:
        if infer_GT == "0|0":
            FN = 1
            TN = 1
        elif infer_GT in ["1|0","0|1"]:
            TP = 1
            TN = 1
        elif infer_GT == "1|1":
            TP = 1
            FP = 1
    elif real_GT == "1|1":
        if infer_GT == "0|0":
            FN = 2
        elif infer_GT in ["1|0","0|1"]:
            TP = 1
            FN = 1
        elif infer_GT == "1|1":
            TP = 2
    
    return [TP,TN,FP,FN]
        

            


def GT_matching():
    #print()
    real_dic = mk_dic(ref)
    infer_dic = mk_dic(infer)
    #print(real_dic)
    #print(infer_dic)
    #TN = 0
    #TP = 0
    #FN = 0
    #FT = 0
    C_table = [0,0,0,0]
    out.write("ID\tTP\tTN\tFP\tFN\n")
    for i in infer_dic:
        # i = HLA_A....
        C_table = [0,0,0,0]
        if i not in real_dic:
            continue
        else:
            print(i)
            for j in infer_dic[i].keys():
                #print(j)
                if j not in real_dic[i]:
                    continue
                else:
                    #print("here")
                    tmp_table = check_GT(real_dic[i][j],infer_dic[i][j])
                    C_table = [a+b for a,b in zip(C_table,tmp_table)]
        new_line = i + "\t" + '\t'.join([str(s) for s in C_table]) + "\n"
        print(new_line)
        out.write(new_line)
        
    print("done")

print("done")
GT_matching()
out.close()

'''
import os,sys

# python 04.check.GT.py input.txt
fileIn = sys.argv[1]
print(fileIn)
out = fileIn.replace(".txt",".after.txt")
out = open(out,'w')
def main():
    #header = "ID\t"
    df = open(fileIn,"r")

    df = df.readlines()
    header = "\t".join([s.strip().split("=")[0] for s in df[0].split()])
    out.write(header+"\n")
    for i in df:
        tmp = i.split()
        new_line = tmp[0]
        for j in tmp[1:]:
            tmp2 = j.split("=")[1]
            if tmp2 == "0|1":
                new_line= new_line + "\t1|0"
            else:
                new_line= new_line + "\t%s"%tmp2
        out.write("%s\n"%new_line)

main()
out.close()




'''
'''
for i in range(1,6):
    if i == 4:
        continue
    else:
        print(str(i))
    print(str(i))
    
a = {'A':{"a":1,"b":2,"c":3},'B':{"aa":11,"bb":22,"cc":33}}
a = {}
"a" in a

for i in a:
    for j in a[i].values():
        print(j)
        print(i)
        print(a[i])
#        print(a[i][a[]])
    

'''