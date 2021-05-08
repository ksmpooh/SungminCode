## UKB data 에서 locus 추출
# input data : chr\tpos\tpvalue 로 구분되어 지는 데이터 범위 지정
# python extract.Locus.py [input] [range] [output]
import os,sys,operator

def fileIn(df):
    out = open(df,"r")
    return [s.replace("\n","").strip() for s in out]

def fileOut(df):
    out = open(df,"w")
    for i in out:
        out.write(i + "\n")
    out.close()


def make_df_sorting_withpvalue(dataIn):
    dic = {}
    for i in range(1,22+1):
        dic[i] = {}
    for i in dataIn:
        chr,pos,pvalue = i.split("\t")
        dic[chr].setdefault(int(pos),int(pvalue))

    for i in range(1,22+1):
        dic[str(i)] = sorted(dic[str(i)],key= operator.itemgetter(1))    
    return dic

#def make_df_sorting()


#def top_rightleft(pos_pvalue,range_length):

def simple_x_distance(df,distance):
    #chr    pos
    for i in df:
        chr,pos = i.split("\t")
        
    




def main():
    datain = sys.argv[1]
    lenth = sys.argv[2]
    dataout = sys.argv[3]
    
    df = fileIn(datain)
    out = open(datout,'w')

    #df = make_df_sorting(df)
    #for i in range(1,22+1):


