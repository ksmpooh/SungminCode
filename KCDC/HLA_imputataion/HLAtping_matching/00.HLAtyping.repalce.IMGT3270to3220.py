'''R 전처리 작업
library(tidyverse)
library(stringr)
library(readxl)
setwd("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/")


df <-read_excel("HLA.typing.Final.result_529sample_IMGT3370.xlsx")
head(df)
ref <-read_excel("HLAtype.check.withID_IMGT3320_genotype 결과.xlsx")


head(ref)


write.table(df,"HLA.typing.Final.result_529sample_IMGT3370.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(ref,"HLAtype.check.withID_IMGT3320_genotype.txt",col.names = T,row.names = F,quote = F,sep = "\t")
'''

# 20230427
# HLA typing data IMGT 3370 to 3320
## df "HLA.typing.Final.result_529sample_IMGT3370.txt"
'''
  Sample ID     NGS_A.1 NGS_A.2 NGS_B.1 NGS_B.2 NGS_C.1 NGS_C.2 NGS_D…¹ NGS_D…² NGS_D…³ NGS_D…⁴ NGS_D…⁵ NGS_D…⁶ NGS_D…⁷ NGS_D…⁸ NGS_D…⁹
  <chr>  <chr>  <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
1 CDC022 NIH19… A*02:0… A*33:0… B*44:0… B*40:0… C*03:0… C*14:0… DQA1*0… DQA1*0… DQB1*0… DQB1*0… DPA1*0… DPA1*0… DPB1*0… DPB1*0… DRB1*0…
2 CDC023 NIH19… A*02:0… A*02:0… B*15:1… B*40:0… C*01:0… C*03:0… DQA1*0… DQA1*0… DQB1*0… DQB1*0… DPA1*0… DPA1*0… DPB1*0… DPB1*0… DRB1*0…
3 CDC024 NIH19… A*02:0… A*33:0… B*40:0… B*51:0… C*03:0… C*03:0… DQA1*0… DQA1*0… DQB1*0… DQB1*0… DPA1*0… DPA1*0… DPB1*0… DPB1*0… DRB1*1…
4 CDC025 NIH19… A*33:0… A*24:0… B*58:0… B*54:0… C*01:0… C*03:0… DQA1*0… DQA1*0… DQB1*0… DQB1*0… DPA1*0… DPA1*0… DPB1*0… DPB1*0… DRB1*0…
5 CDC026 NIH19… A*02:0… A*02:0… B*38:0… B*27:0… 
'''
###ref "HLAtype.check.withID_IMGT3320_genotype.txt"
'''
    0      1    2           3           4        
HLAID	gene	IMGT3370	IMGT3320	IMGT3320 분석 결과_추가
2020KDCA001	A	A*02:01:01:50	NA	A*02:01:01:01
2020KDCA001	C	C*07:02:01:32	NA	C*07:02:01:03
2020KDCA003	DQA1	DQA1*05:05:01:15	NA	DQA1*05:05:01:03
2020KDCA004	DQA1	DQA1*05:05:01:15	NA	DQA1*05:05:01:03
2020KDCA005	B	B*15:02:01:01	B*15:02:01	B*15:02:01
2020KDCA006	B	B*35:01:01:13	NA	B*35:01:01:05
2020KDCA008	B	B*44:03:01:13	NA	B*44:03:01:01
2020KDCA008	C	C*14:03:01	C*14:03	C*14:03
2020KDCA009	DPA1	DPA1*02:01:01:06	NA	DPA1*02:01:01:01
'''


import os,glob


wDir = "/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/"

def main():
    out = open(wDir + "HLA.typing.Final.result_529sample_IMGT3320convert.txt","w")
    
    ref = open(wDir + "HLAtype.check.withID_IMGT3320_genotype.txt","r")
    ref = ref.readlines()
    ref.pop(0)

    df = open(wDir + "HLA.typing.Final.result_529sample_IMGT3370.txt","r")

    df = df.readlines()
    out.write(df.pop(0))

    df_dic = {}
    print("mk dic...")
    for i in df:
        tmp = i.split("\t")
        df_dic[tmp[0]] = i.strip()
    
    print("check and replace...")
    for i in ref:
        tmp = i.split("\t")
        HLAID = tmp[0]
        IMGT3370 = tmp[2].strip()
        IMGT3320 = tmp[4].strip()
        if HLAID in df_dic.keys():
            if IMGT3320 == "-":
                IMGT3320 = "0"
            print("HLA ID ( %s ) : %s -> %s"%(HLAID,IMGT3370,IMGT3320))
            df_dic[HLAID] = df_dic[HLAID].replace(IMGT3370,IMGT3320)
            #print(df_dic[HLAID])
        else:
            continue
    
    for i in df_dic:
        out.write("%s\n"%(df_dic[i]))

    out.close()

    print("done....")

main()


    




