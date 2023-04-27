#### 20230323 20230323 KOTRY 회의

####/BDATA/smkim/KR_allogenomics/04.immune_cell_signal/

# co-signal_pairtable.txt
#0          1       2           3           4        5                      6       7   
#Celltype	Signal1	Co-signals	Receptor	Ligand	Intracellular_signaling	RLpair	RLIpair
#T cell (αβ)	TCR	Costimulator	CD28	CD80,CD86,ICOSLG	PIK3CD, GRB2,  VAV1	CD28,CD80,CD86,ICOSLG	CD28,CD80,CD86,ICOSLG,PIK3CD, GRB2,  VAV1

'''KR.KD.immune.cell.co-signal_targetGene.alleleMatching01.Score_Sum.txt
KBA_ID.KD	KBA_ID.KR	AKT1_Score_SUM	BIRC5_Score_SUM	BTLA_Score_SUM
NIH19KT0267	NIH19KT0023	0	1	0
NIH19KT0268	NIH19KT0024	0	0	0
NIH19KT0269	NIH19KT0025	0	0	2
NIH19KT0270	NIH19KT0026	0	1	0
'''
'''
co-signal_pairtable_withINDEX.txt
0       1       2               3        4
Signal1	RLpair	RLpair_index	RLIpair	RLIpair_index
TCR	CD28,CD80,CD86,ICOSLG	RLpair_1	CD28,CD80,CD86,ICOSLG,PIK3CD, GRB2,  VAV1	RLIpair_1
TCR	ICOS,ICOSLG	RLpair_2	ICOS,ICOSLG,PIK3CD	RLIpair_2
TCR	TREML2,CD276	RLpair_3	TREML2,CD276,NA	RLIpair_3
TCR	NA,CD274,PDCD1LG2	RLpair_4	NA,CD274,PDCD1LG2,NA	RLIpair_4
'''
import os,glob
wDir = "/BDATA/smkim/KR_allogenomics/04.immune_cell_signal/"
def main():
    datain = open(wDir + "KR.KD.immune.cell.co-signal_targetGene.alleleMatching01.Score_Sum.txt","r")
    refin = open(wDir + "co-signal_pairtable_withINDEX.txt","r")    
    outdata = open(wDir + "KR.KD.immune.cell.co-signal_targetGene_RLpair.alleleMatching01.Score_Sum.txt","w")
    #out.write("KBA_ID.KD\tKB_ID.KR\t")
    #out = [["KBA_ID.KD","KB_ID.KR"]]
    out =[]
    df = datain.readlines()
    header = df[0].replace("\n","").replace("_Score_SUM","").split("\t")
    print(header)
    df = [i.strip().split("\t") for i in df]
    for i in df:
        id1 = i[0]
        id2 = i[1]
        tmp = [id1,id2]
        out.append(tmp)
        
    ref_header = refin.readline()
    while 1:
    #for i in range(0,10):
        line = refin.readline().replace("\n","")
        if not line:
            break
        tmp = line.split("\t")
        RLpair = tmp[1].split(",")
        RLpair_index = tmp[2]

        df_idx = []

        for gene in RLpair:
            gene = gene.strip()
            if (gene == "NA") or (gene not in header):
                continue
            df_idx.append(header.index(gene))
        #print(df_idx,RLpair)
        out[0].append(RLpair_index)
        for r_idx in range(1,len(out)):
            scores = []
            for score_idx in df_idx:
                scores.append(int(df[r_idx][score_idx]))
            out[r_idx].append(str(sum(scores)))
        print(RLpair_index,RLpair,df_idx,scores)
        #print(RLpair)
        #print(scores)


        RLIpair = tmp[3].split(",")
        RLIpair_index = tmp[4]
        df_idx = []
        for gene in RLIpair:
            gene = gene.strip()
            if (gene == "NA") or (gene not in header):
                continue
            df_idx.append(header.index(gene))
        #print(df_idx,RLpair)
        out[0].append(RLIpair_index)
        for r_idx in range(1,len(out)):
            scores = []
            for score_idx in df_idx:
                scores.append(int(df[r_idx][score_idx]))
            out[r_idx].append(str(sum(scores)))
        print(RLIpair_index,RLIpair,df_idx,scores)
        #print(RLIpair)
        #print(scores)



        #print(RLpair)
        #RLIpair = tmp[3].split(",")
        #RLIpair_index = tmp[4]
    #print(out[0:20])
    for i in out:
        outdata.write("\t".join(i))
        outdata.write("\n")
    outdata.close()


main()






    

