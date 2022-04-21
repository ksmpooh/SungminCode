###1 make data frame from .raw 
###2 pair table score
## python3 .py [score table input.txt]
## python3 .py [score table input.txt]
## python3 .py [_raw]
import os,glob,sys,re

inData = sys.argv[1]
#outData = inData.replace(".txt","_dataproScore.txt")

##Rscript --vanilla 03.ploting.R [input.raw] [output.txt]
wDir = os.getcwd()
def give_score(a):
    if a.count("NA") > 0:
        return "\t0"
    #print(a)
    if a[0] == a[1]:
        return "\t2"
    elif abs(int(a[0]) - int(a[1])) == 1:
        return "\t1"
    else:
        return "\t0"


def make_score_table(pair_table):
    print("Make_pair table....from .raw file")
    print("FILE IN : "+pair_table)
    outData = pair_table.replace(".raw","_pairTable.txt")
    os.system("Rscript --vanilla 03.data.Pair.matching.R %s %s"%(inData,outData))
    return outData


def calculate_score(dataIn,theme):
    print("Calculate Score : ..."+dataIn)
    print("theme : "+theme)
    df = open(dataIn,"r")
    outData = dataIn.replace("_pairTable","_ScoreTable")
    out = open(outData,"w")
    ori_header= df.readline().replace("\n","").split("\t")
    print(ori_header[0:5])
    new_header = "KBA_ID.KD\tKBA_ID.KR"


    refs = []
    for ref in ori_header[2:]:
        if ref[-1] == 'y':
            break
        refs.append(ref.replace(".x",""))
        #new_header = new_header + "\t" + ref[0].replace("X","").replace(".x","")
        new_header = new_header + "\t" + ref.replace(".x","")
    #new_header = new_header + "\n"
    new_header = new_header + "\t%s_Score_SUM\n"%theme
    out.write(new_header)

    ref_dict = {}
    print("ref len : "+str(len(refs)))
    for ref in refs:
        a = [x for x in ori_header if ref in x]
        ref_dict.setdefault(ref,[ori_header.index(a[0]),ori_header.index(a[1])])

    while 1:
        line = df.readline()
        if not line:
            break
        line =line.replace("\n","").split("\t")
        out.write("%s\t%s"%(line[0],line[1]))
        score = ""
        for ref in refs:
            i1,i2 = ref_dict[ref]
            #print(i1,i2)
            out.write(give_score([line[i1],line[i2]]))
            score = score + give_score([line[i1],line[i2]])
        score_Sum = sum([int(s) for s in score[1:].split("\t")])
        #out.write("\n")
        out.write("\t%s\n"%(str(score_Sum)))
    print("Done..calculate_score")
    out.close()



def main():
    print("main......")
    #KR.KD.transmembrane.raw
    theme = inData.replace(wDir +"/","").replace("KR.KD.","").replace(".raw","")

    tmp_df = make_score_table(inData)
    calculate_score(tmp_df,theme)
    print("Done! ALL python script...")


main()


    
