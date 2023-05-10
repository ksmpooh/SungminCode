## HLa eplet matching data inout
## HLA input
##
## python 00.eplet.matching.excel.inout.py [class 1 or 2] [input] 
## ex) python 00.eplet.matching.excel.inout.py 1 HLA_classI for eplet
## class I/II
## input : no header file 
import os,sys,glob
import pyxlsb
import openpyxl as oxl
import pandas as pd


#mhc_class = sys.argv[1]
#hla_df = sys.argv[2]

mhc_class = "1"
hla_df = "/Users/ksmpooh/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/KR.KD.HLAimp.foreplet.class1.xlsx"
#### INPUT
'''
#classI II  file example)
RecInfo	Rec1stA	Rec2ndA	Rec1stB	Rec2ndB	Rec1stC	Rec2ndC	DonorInfo	Don1stA	Don2ndA	Don1stB	Don2ndB	Don1stC	Don2ndC
ID	1stDRB	2ndDRB	1stDRW	2ndDRW	1stDQB	2ndDQB	1stDQA	2ndDQA	1stDPB	2ndDPB	1stDPA	2ndDPA	
ID	1stDRB	2ndDRB	1stDRW	2ndDRW	1stDQB	2ndDQB	1stDQA	2ndDQA	1stDPB	2ndDPB	1stDPA	2ndDPA
'''

'''class I
NIH19KT6603	A*24:02	A*29:01	B*07:05	B*67:01	C*07:02	C*15:05	NIH19KT6653	A*24:02	A*02:06	B*39:01	B*55:02	C*03:04	C*01:02
NIH19KT6606	A*11:01	A*11:01	B*46:01	B*46:01	C*01:02	C*01:02	NIH19KT6654	A*11:01	A*11:01	B*46:01	B*46:01	C*01:02	C*01:02
'''
''' class II
NIH19KT6603	DRB1*04:10	DRB1*08:03	x	x	DQB1*04:02	DQB1*03:01	DQA1*03:03	DQA1*06:01	DPB1*05:01	DPB1*05:01	DPA1*02:02	DPA1*02:02	NIH19KT6653	DRB1*04:05	DRB1*14:54	x	x	DQB1*05:02	DQB1*04:01	DQA1*01:04	DQA1*03:03	DPB1*04:02	DPB1*04:02	DPA1*01:03	DPA1*01:03
NIH19KT6606	DRB1*08:03	DRB1*08:03	x	x	DQB1*06:01	DQB1*06:01	DQA1*01:03	DQA1*01:03	DPB1*02:01	DPB1*02:02	DPA1*01:03	DPA1*02:02	NIH19KT6654	DRB1*08:03	DRB1*08:03	x	x	DQB1*06:01	DQB1*06:01	DQA1*01:03	DQA1*01:03	DPB1*02:02	DPB1*02:02	DPA1*02:02	DPA1*02:02
'''


def main():
    print("INPUT MHC class II : %s"%mhc_class)
    if mhc_class not in ["1","2"]:
        print("MCH class only for 1 or 2")
        print("Please INPUT correctly (1 or 2)")
    else:
        print("INPUT MHC class II : %s"%mhc_class)
        if mhc_class == "1":
            hlamatchmaker = "/Users/ksmpooh/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/20230420_eplet/ABC_Eplet_Matching_4.0.xlsb"
        else:
            hlamatchmaker = "/Users/ksmpooh/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/20230420_eplet/DRDQDP_Eplet_Matching_3.1.xlsb"

    '''
    load_df = oxl.load_workbook("%s"%hla_df)
    df_pre = load_df['Sheet1']
    #print(df_pre[0:2])
    df = []
    for row in df_pre:
        tmp = []
        for cell in row:
            tmp.append(cell)
        df.append(tmp)
    '''
    df = pd.read_excel(hla_df,engine='openpyxl',header=None)
    
    print(df[0:5])
    print(df[0,1])

main()

    

    





