### infer vs NGS
## DATA check

## python python.py [real] [infer] [digit] [output]
## digit : 2,4
'''infer
ID	HLA-A_1	HLA-A_2	HLA-B_1	HLA-B_2	HLA-C_1	HLA-C_2	HLA-DRB1_1	HLA-DRB1_2	HLA-DPA1_1	HLA-DPA1_2	HLA-DPB1_1	HLA-DPB1_2	HLA-DQA1_1	HLA-DQA1_2	HLA-DQB1_1	HLA-DQB1_2
NIH19KT2304	33:03:23	24:02:01:03	51:01:05	58:01:07	14:02:01	03:02:02:01	11:06:01	14:10	01:03:01:02	02:02:02	05:01:01	04:02:01:02	01:04:01:02	05:05:01:02	03:01:01:03	05:02:01
'''

'''NGS
'''


import os,sys,glob
'''
b={}
b["ID1"] = {}
b["ID1"]["A"] = []
b["ID1"]["A"].append(1)
b["ID1"]["A"].append(2)
print(b)
'''

def digit_split(hlatype,digit):
    #print(hlatype)
    if hlatype == "NA":
        return "NA"
    else:
        tmp = hlatype.split(":")
        if digit == "2":
            return tmp[0]
        elif digit == "4":
            return tmp[0] + ":" + tmp[1]
        else:
            return "NA"


'''
{ID1:{HLA-A:[01,02],HLA-B:[02,08]},
 ID2:{HLA-A:[01,03],HLA-B:[02,04]},
 ...
 }
'''
def mk_dic_forMatching(df,genes,digit):
    print("mk_dic...")
    new_dic = {}
    header = df.pop(0)
    header = header.split()
    #print(header)
    for i in df:
        #print(i)
        i = i.strip().split()
        #print(i)
        new_dic[i[0]] = {}
        for gene in genes:
            new_dic[i[0]][gene] = []
            #print(i[header.index[gene+"_1"]])
            #new_dic[i[0]][gene].append(i[header.index[gene+"_1"]])
            new_dic[i[0]][gene].append(digit_split(i[header.index(gene+"_1")],digit))
            #new_dic[i[0]][gene].append(i[header.index[gene+"_2"]])
            new_dic[i[0]][gene].append(digit_split(i[header.index(gene+"_2")],digit))
    return new_dic


'''
    header = ["ID","A.match","A.wrong","A.empty","B.match","B.wrong","B.empty","C.match","C.wrong","C.empty",
           "DRB1.match","DRB1.wrong","DRB1.empty","DPA1.match","DPA1.wrong","DPA1.empty",
           "DPB1.match","DPB1.wrong","DPB1.empty","DQA1.match","DQA1.wrong",
           "DQA1.empty","DQB1.match","DQB1.wrong","DQB1.empty","match_Sum","wrong_Sum","empty_Sum"]
'''

def matching_table(real,infer,genes):
    print("matching_table.....")
    header = ["ID"]
    for gene in genes:
        for i in [gene+".match",gene+".wrong",gene+".empty"]:
            header.append(i)
    header.append("match_Sum")
    header.append("wrong_Sum")
    header.append("empty_Sum")

    out = []
    out.append(header)
    for ID in infer:
        tmp = []
        tmp.append(ID)
        real_ID = real[ID]
        match_sum = 0
        wrong_sum = 0
        empty_sum = 0
        for gene in genes:
            match = 0
            wrong = 0
            empty = 0
            if (infer[ID][gene].count("NA") == 2) | (real[ID][gene].count("NA") == 2):
                empty = 2
            elif (infer[ID][gene].count("NA") == 1) | (real[ID][gene].count("NA") == 1):
                empty = 1
            #life = empty

            for i in infer[ID][gene]:
                if match + wrong + empty == 2:
                    #print("life : %s %s %s"%(match,wrong,empty))
                    break
                elif i == "NA":
                    continue
                else:
                    if i in real_ID[gene]:
                        match = match + 1
                        real_ID[gene].remove(i)
                    #elif (i not in ngs) & (life ==0):
                    else: #
                        wrong = wrong + 1
            match_sum = match_sum + match
            wrong_sum = wrong_sum + wrong
            empty_sum = empty_sum + empty
        
            tmp.append(str(match))
            tmp.append(str(wrong))
            tmp.append(str(empty))
        tmp.append(str(match_sum))
        tmp.append(str(wrong_sum))
        tmp.append(str(empty_sum))
        print("print tmp")
        print(tmp)
        out.append(tmp)
    return out
    





    

'''
def match_process(ori,gene):
    df = ori
    
    #print(df)
    for i in [gene+".match",gene+".wrong",gene+".empty"]:
        df[0].append(i)
    #print(df) 

    for index,i in enumerate(df[1:]):
        #print(index,i)
        match = 0
        wrong = 0
        empty = 0
        sm = i[0:1+1]
        ngs = i[2:3+1]
        #print(ngs)
        if (sm.count("NA") == 1) | (ngs.count("NA") == 1):
            empty = 1
        elif (sm.count("NA") == 2) | (ngs.count("NA") ==2):
            empty = 2 
        #print("count : empty %s"%str(empty))
        life = empty
        for a in sm:
            if match + wrong + empty == 2:
            #    print("Break")
                break
            if a in ngs:
            #    print("match")
                match = match + 1
                ngs.remove(a)
            elif (a not in ngs) & (life == 0):
            #    print("wrong")
                wrong = wrong +  1
            elif (a not in ngs) & (life == 1):
                life = life -  1
        #print("last")
        df[index+1].append(match)
        df[index+1].append(wrong)
        df[index+1].append(empty)
        #print(df[index+1])
    return df

def gene_df(ori,a1,a2,a3,a4):
    outdf = []
    for i in ori:
        x = []
        x.append(i[int(a1)])
        x.append(i[int(a2)])
        x.append(i[int(a3)])
        x.append(i[int(a4)])
        outdf.append(x)
        #print(x)
    return outdf 

def final_processing(datain,dataout):
    df = filein(datain)
    df = [s.split("\t") for s in df]
    df = match_process(df)
    fileout(df,dataout)
'''


def main():
    real_fileName = sys.argv[1]
    infer_fileName = sys.argv[2]
    digit = sys.argv[3]
    fileOut = sys.argv[4]

    real_df = open(real_fileName,"r")
    infer_df = open(infer_fileName,"r")

    real_df = real_df.readlines()
    infer_df = infer_df.readlines()

    genes = ["HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
    genes = genes

    read_dic = mk_dic_forMatching(real_df,genes,digit)
    infer_dic = mk_dic_forMatching(infer_df,genes,digit)

    print(read_dic)
    print(infer_dic)

    out = matching_table(read_dic,infer_dic,genes)
    print(out)

    with open(fileOut,'w') as out_write:
        for i in out:
            out_write.write("%s\n"%("\t".join(i)))
    print("done")
    

main()

#python test.py merge.df.txt HLA.scan.NIH19KT1012.result_processing.txt 2 test.txt
