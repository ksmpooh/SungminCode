### infer vs NGS
## DATA check

## python python.py [real] [infer] [digit] [output]
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

def digit_split(df,digit):
    if digit not in ["2","4","6","8"]:
        print("Digit is not in 2,4,6,8")
        return
    elif digit == "2":
        a= []
    elif digit == "4":
        a = []
    elif digit == "6":
        a =[]
    else:
        a =[]
    return a


def mk_dic_forMatching(df,genes,digit):
    new_dic = {}
    header = df.pop(0)
    for i in df:
        i = i.strip().split()
        new_dic[i[0]] = {}
        for gene in genes:
            new_dic[i[0]][gene] = []
            #new_dic[i[0]][gene].append(i[header.index[gene+"_1"]])
            new_dic[i[0]][gene].append(digit_split(i[header.index[gene+"_1"]],digit))
            #new_dic[i[0]][gene].append(i[header.index[gene+"_2"]])
            new_dic[i[0]][gene].append(digit_split(i[header.index[gene+"_2"]],digit))
    return new_dic



def matching_table(real,infer):
    print(i)
    




def main():
    real_fileName = sys.argv[1]
    infer_fileName = sys.argv[2]
    digit = sys.argv[3]

    real_df = open(real_fileName,"r")
    infer_df = open(infer_fileName,"r")

    real_df = real_df.readlines()
    infer_df = infer_df.readlines()

    genes = ["HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
    genes = genes

    read_dic = mk_dic_forMatching(real_df,genes)
    infer_dic = mk_dic_forMatching(infer_df,genes)
