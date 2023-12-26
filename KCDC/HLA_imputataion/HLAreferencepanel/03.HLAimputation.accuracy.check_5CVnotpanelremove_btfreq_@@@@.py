### infer vs NGS
## DATA check

## python python.py [real] [infer] [digit] [outDir]
## digit : 0 : 자동
##         2 : td
#          4 : fd     
#          42 : fdtotd
'''
1) michigan HLA imputation result
ID      HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2
NIH19KT0019     A*33    A*24    B*07    B*58
NIH19KT0013     A*24    A*02    B*40    B*15
NIH19KT0008     A*03    A*26    B*44    B*15
NIH19KT0254     A*02    A*33    B*07    B*46
NIH19KT0264     A*03    A*02    B*44    B*37
NIH19KT0410     A*02    A*68    B*27    B*38

ID      HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2
NIH19KT0019     A*33:03:01:01   A*24:02:01:01   B*07:02:01:01   B*58:01:01:01
NIH19KT0013     A*24:02:01:01   A*24:86N        B*40:02:01:01   B*15:01:01:01
NIH19KT0008     A*03:01:01:01   A*26:01:01:01   B*44:02:01:01   B*15:01:01:01
NIH19KT0254     A*33:03:01:01   A*02:07:01      B*07:02:01:01   B*46:01:01

'''

'''NGS

1) Engenbio HLA accutest KIT
/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtyping.vs.HLAimp.accuracy.check.NGStyping.DATA.txt
      1 FID     IID     pID     mID     SEX     PHENO   HLA_A.1 HLA_A.2
      2 NIH19KT0406     NIH19KT0406     0       0       0       0       A*02:01:149     A*33:03:01:01
      3 NIH19KT0407     NIH19KT0407     0       0       0       0       A*02:01:01:01   A*02:01:01:01
      4 NIH19KT0408     NIH19KT0408     0       0       0       0       A*02:01:01:01   A*33:03:01:01
      5 NIH19KT0409     NIH19KT0409     0       0       0       0       A*33:03:01:01   A*24:02:01:01
      6 NIH19KT0410     NIH19KT0410     0       0       0       0       A*02:06:01:01   A*02:03:01
      7 NIH19KT0411     NIH19KT0411     0       0       0       0       A*02:06:01:01   A*31:01:02:01
      8 NIH19KT0412     NIH19KT0412     0       0       0       0       A*02:01:01:01   A*02:01:01:01
      9 NIH19KT0413     NIH19KT0413     0       0       0       0       A*02:01:149     A*33:03:01:01
     10 NIH19KT0414     NIH19KT0414     0       0       0       0       A*11:01:01:01   A*26:01:01:01
     11 NIH19KT0415     NIH19KT0415     0       0       0       0       A*02:01:01:01   A*02:01:01:01
     12 NIH19KT0416     NIH19KT0416     0       0       0       0       A*24:02:01:01   A*30:01:01:01
     13 NIH19KT0973     NIH19KT0973     0       0       0       0       A*24:02:01:01   A*24:02:01:01
     14 NIH19KT0974     NIH19KT0974     0       0       0       0       A*33:03:01:01   A*24:02:01:01
     15 NIH19KT0975     NIH19KT0975     0       0       0       0       A*11:01:01:01   A*31:01:02:01
     16 NIH19KT0976     NIH19KT0976     0       0       0       0       A*02:01:01:01   A*11:01:01:01

2) HLA-TAPAS nomen cleaner -> header 추가 후 사용!!!
/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLA.4digit.520sample.nomenclean.chped
      1 NIH19KT0406     NIH19KT0406     0       0       0       0       A*02:01 A*33:03
      2 NIH19KT0407     NIH19KT0407     0       0       0       0       A*02:01 A*02:01
      3 NIH19KT0408     NIH19KT0408     0       0       0       0       A*02:01 A*33:03
      4 NIH19KT0409     NIH19KT0409     0       0       0       0       A*33:03 A*24:02
      5 NIH19KT0410     NIH19KT0410     0       0       0       0       A*02:06 A*02:03
      6 NIH19KT0411     NIH19KT0411     0       0       0       0       A*02:06 A*31:01
      7 NIH19KT0412     NIH19KT0412     0       0       0       0       A*02:01 A*02:01
      8 NIH19KT0413     NIH19KT0413     0       0       0       0       A*02:01 A*33:03
      9 NIH19KT0414     NIH19KT0414     0       0       0       0       A*11:01 A*26:01
     10 NIH19KT0415     NIH19KT0415     0       0       0       0       A*02:01 A*02:01
     11 NIH19KT0416     NIH19KT0416     0       0       0       0       A*24:02 A*30:01
     12 NIH19KT0973     NIH19KT0973     0       0       0       0       A*24:02 A*24:02
     13 NIH19KT0974     NIH19KT0974     0       0       0       0       A*33:03 A*24:02
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
real_fileName = sys.argv[1]
print("Real_file : %s"%real_fileName)
if ".chped" in real_fileName:
    real_type = "Nomencleaner"
else:
    real_type = "RealNGStyping"


infer_fileName = sys.argv[2]
print("HLAimp_file : %s"%infer_fileName)
digit = sys.argv[3]

outDir = sys.argv[4]
fileOut = outDir + "/" + infer_fileName.split("/")[-1].replace(".txt",".cmp_%s.txt"%real_type)

if digit not in ["0","2","4","42"]:
    print("Digit Error (0, 2, 4, 42)")
    sys.exit("python python.py [real] [infer] [digit] [output]")
elif digit == "0":
    if "_td.txt" in infer_fileName:
        print("Digit 0 : file name with td -> digit = 2")
        digit = '2'
    else:
        print("Digit 0 : file name with fd -> digit = 4")
        digit = '4'
elif digit == "42":
    fileOut = fileOut.replace(".txt",".fdvstd.txt")
    digit = "2"
    print("Digit %s"%(digit))
else:
    print("Digit %s"%(digit))



#fileOut = sys.argv[4]
miss_matching_fileName = fileOut.replace(".txt","_missINFO.txt")


'''
def digit_split(hlatype,digit):
    #print(hlatype)
    if hlatype == "NA" or hlatype == "0":
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
### not panel remove

def digit_split(hlatype,digit,panel_type):
    #print(hlatype)
    if hlatype == "NA" or hlatype == "0":
        return "NA"
    else:
        tmp = hlatype.split(":")
        if digit == "2":
            if tmp[0] not in panel_type:
                return "NA"
            else:
                return tmp[0]
        elif digit == "4":
            if tmp[0] + ":" + tmp[1] not in panel_type:
                return "NA"
            else:
                return tmp[0] + ":" + tmp[1]
        else:
            return "NA"
###


'''
{ID1:{HLA-A:[01,02],HLA-B:[02,08]},
 ID2:{HLA-A:[01,03],HLA-B:[02,04]},
 ...
 }
'''
def mk_dic_forMatching(df,genes,digit,panel_type):
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
            new_dic[i[0]][gene].append(digit_split(i[header.index(gene+".1")],digit,panel_type))
            #new_dic[i[0]][gene].append(i[header.index[gene+"_2"]])
            new_dic[i[0]][gene].append(digit_split(i[header.index(gene+".2")],digit,panel_type))
    return new_dic


'''
    header = ["ID","A.match","A.wrong","A.empty","B.match","B.wrong","B.empty","C.match","C.wrong","C.empty",
           "DRB1.match","DRB1.wrong","DRB1.empty","DPA1.match","DPA1.wrong","DPA1.empty",
           "DPB1.match","DPB1.wrong","DPB1.empty","DQA1.match","DQA1.wrong",
           "DQA1.empty","DQB1.match","DQB1.wrong","DQB1.empty","match_Sum","wrong_Sum","empty_Sum"]
'''

def matching_table(real,infer,genes):
    mm = open(miss_matching_fileName,'w')
    mm.write("ID\tgene\tReal\tImp\n")

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
                        mm.write("%s\t%s\t%s\t%s\n"%(ID,gene,real_ID[gene][0],i))
            match_sum = match_sum + match
            wrong_sum = wrong_sum + wrong
            empty_sum = empty_sum + empty
        
            tmp.append(str(match))
            tmp.append(str(wrong))
            tmp.append(str(empty))
        tmp.append(str(match_sum))
        tmp.append(str(wrong_sum))
        tmp.append(str(empty_sum))
        #print("print tmp")
        #print(tmp)
        out.append(tmp)
    mm.close()
    return out

'''
Gene	value	n	prop	freq
A	A*01:01	19	0.0182692307692308	less_common
A	A*02:01	174	0.167307692307692	common
A	A*02:03	6	0.00576923076923077	rare
A	A*02:05	1	0.000961538461538462	rare
A	A*02:06	94	0.0903846153846154	common
'''

def main():
#    real_fileName = sys.argv[1]
#    infer_fileName = sys.argv[2]
#    digit = sys.argv[3]
#    fileOut = sys.argv[4]

    freq_ref = open("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/KMHC_HLAtype_frequency.txt","r")
    freq_ref = freq_ref.readlines()
    common = []
    less_common = []
    rare = []
    for i in freq_ref[1:]:
        tmp = i.split()
        if tmp[4] == "common":
            common.append(tmp[1])
        elif tmp[4] == "rare":
            rare.append(tmp[1])
        else:
            less_common.append(tmp[1])

        


    real_df = open(real_fileName,"r")
    infer_df = open(infer_fileName,"r")

    real_df = real_df.readlines()
    infer_df = infer_df.readlines()

    ##### real_df에 없는 type 지우기 위한 체크리스트
    infer_sample_id = []
    for i in infer_df[1:]:
        tmp = i.split()
        infer_sample_id.append(tmp[0])

    #real_df_tmp = real_df
    #real_df = []
    #real_df.append(real_df_tmp.pop(0))
    panel_type = []
    target_type = []

    for i in real_df[1:]:
        tmp = i.strip().split()
        if tmp[0] in infer_sample_id: 
            ## imp에만 있는 type 추출
            for j in tmp[2:]:
                if j in ["0","NA"]:
                    pass
                else:
                    target_type.append(j)
        else:
            ## panel에만 있는 type 추출
            for j in tmp[2:]:
                if j in ["0","NA"]:
                    pass
                else:
                    panel_type.append(j)

    panel_type = list(set(panel_type)) # 중복제거
    target_type = list(set(target_type)) # 중복제거

    panel_type = [s for s in panel_type if s in target_type]


    #print(panel_type)

    if digit == "2":
        panel_type = [a.split(":")[0] for a in panel_type]
        panel_type = list(set(panel_type)) # 중복제거
    else:
        panel_type = [a.split(":")[0] + ":" + a.split(":")[1] for a in panel_type]
        panel_type = list(set(panel_type)) # 중복제거

    print(panel_type)
    print(len(panel_type))
    #####

    genes = ["HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1"]
    genes = genes

    header = infer_df[0].strip().split()
    new_genes = []
    for gene in genes:
        if gene +".1" in header:
            new_genes.append(gene)
        else:
            pass
    genes = new_genes

    read_dic = mk_dic_forMatching(real_df,genes,digit,panel_type)
    print(read_dic)
    infer_dic = mk_dic_forMatching(infer_df,genes,digit,panel_type)
    print(infer_dic)

    #print(read_dic)
    #print(infer_dic)

    out = matching_table(read_dic,infer_dic,genes)
    #print(out)

    with open(fileOut,'w') as out_write:
        for i in out:
            out_write.write("%s\n"%("\t".join(i)))
    print("done")

    

main()

#python test.py merge.df.txt HLA.scan.NIH19KT1012.result_processing.txt 2 test.txt




'''
datain = open("HLA.type.forHLATAPA.txt","r")
dataout = open("HLA.type.forHLATAPA_2field.txt","w")

header = datain.readline()
dataout.write(header)
while 1:
    line = datain.readline()
    if not line:
        break
    tmp = line.split()
    for i in range(6,len(tmp)):
        tmp[i] = ":".join(tmp[i].split(":")[0:2])
    dataout.write("\t".join(tmp))

dataout.close()
'''