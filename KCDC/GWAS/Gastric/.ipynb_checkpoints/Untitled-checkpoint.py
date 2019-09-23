#!/usr/bin/env python
# coding: utf-8

# ### ref file match and change 
# plink bim 파일에서 A1, A2를 ref 파일내에 있는 allele과 비교하여 매칭되는 allele이 ref allele입니다.
# 매칭 후 chr:pos_ref/alt 로 bim 파일내에 있는 affy id 를 변환하여
# 1KG P3와 공통인 마커만 추려 합친 후 PCA 분석을 진행하면 됩니다.
# <분석과정>
# 1. plink bim 파일 A1, A2에서 ref allele 확인
# 2. bim 파일 내 affy id를 new id(chr:pos_ref/alt)로 변환
# 3. 1KG P3와 공통인 마커 추린 후 두 데이터를 합치기
# 4. PCA 분석하기
# 5. 인종별 색을 다르게 하여 plot 그리기 (1000GP_Phase3.sample에서 인종별 ID 확인)



import os
import pandas as pd




wdir = "c:/Users/user/Desktop/KCDC/Gastric/Ref/"



def fileRead(fileIn):
    f = open(fileIn,'r')
    inData = [r.replac("r","").replace("\n","").replace("\t","") for r in f]
    return inData

bim = pd.read_csv(wdir+"Case_Control_merge_rmfreq.bim",delim_whitespace = True,header = None)
ref = pd.read_csv(wdir+"Axiom_KOR.annot.extract.addINDEL.Final.REF.txt",delim_whitespace = True,header = None)



new_id = []
for i in range(0,len(bim)):
    index = ref.loc[ref[0] == bim[1][i]]
    if (index[1] == bim[4][i]).bool():
        new_id.append(str(bim[0][i]) + ":" + str(bim[3][i]) +"_"+str(bim[4][i])+"/"+str(bim[5][i]))
    else:
        new_id.append(str(bim[0][i]) + ":" + str(bim[3][i]) +"_"+str(bim[5][i])+"/"+str(bim[4][i]))
	tmp = bim[4][i]
	bim[4][i] = bim[5][i]
	bim[5][i] = tmp
bim.loc[:,1] = new_id
bim.to_csv(wdir+"match_ref_merge.bim",header = False,index = False, sep = '\t')






