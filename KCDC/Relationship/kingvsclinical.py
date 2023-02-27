import os, sys

# python input.py king.txt clinicarl.txt

'''ref
ID1 ID2 type
NID1    NID2    FS
NID3    NID4    PS

JG.not.predicted.prunedSNP.kin0
0       1   2       3   4       5       6       7       8       9       10      11      12      13
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType
NIH19KT0010	NIH19KT0010	NIH19KT6345	NIH19KT6345	111279	0.0836	0.0162	0.2767	0.2028	0.1315	0.4729	0.0000	0.2365  2nd


'''                                                     


def main():
    df = open(sys.argv[1],"r")
    ref = open(sys.argv[2],"r")
    

    ref = [s.replace("\n","").split() for s in ref]
    #{ID1 : [ID2,relative]}
    ref_dic = {s[0]:[s[1],s[2]] for s in ref}
    for s in ref:
        ref_dic[s[1]] =  [s[0],s[2]]
    

    out = open(sys.argv[1].replace(".kin0","_reuslt.txt"),"w")


    header = df.readline().replace("\n","")
    df = [s.replace("\n","") for s in df]
    
    out.write(header + "\tRelaType\n")

    for i in df:
        tmp = i.split()
        if tmp[0] in ref_dic.keys():
            if tmp[2] == ref_dic[tmp[0]][0]:
                out.write("%s\t%s\n"%(i,ref_dic[tmp[0]][1]))
            else:
                continue
        else:
            continue
    out.close()


#main()










    


