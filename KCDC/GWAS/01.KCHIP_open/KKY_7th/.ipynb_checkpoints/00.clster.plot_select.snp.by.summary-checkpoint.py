#AFFX-KIT-000001-A

# sumamry : /DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19/AxiomGT1.summary.txt
# ref : /DATA/smkim/KKY/02.1stQC/OUTPUTs/all/clusterplot/5e-8_marker.list.txt
# summary withtout # : /DATA/smkim/KKY/02.1stQC/OUTPUTs/all/clusterplot/summary.without.shap.txt
def main():
    # marker list 안에 2nd column이 마커 ID
    refin = open('/DATA/smkim/KKY/02.1stQC/OUTPUTs/all/clusterplot/5e-8_marker.list.txt','r')
    ref = []
    out = open("/DATA/smkim/KKY/02.1stQC/OUTPUTs/all/clusterplot/summary_for.cluster.txt",'w')
    while 1:
        a = refin.readline()
        if not a: 
            break
        a = a.split()
        ref.append(a[1])
    print(ref)
    refin.close()
    count = 1
    df = open("/DATA/smkim/KKY/01.GenotypeCalling/OUTPUTs/1st.apt-1.19/AxiomGT1.summary.txt","r")
    while True:
        a = df.readline()
        if not a:
            break
        if a[0] == "#": # 주석제거
            continue
        if count == 1: # header 추가
            out.write(a)
            count = 2
        b = a.split() 
        if b[0][:-2] in ref: # marker 비교 ID-A, ID-B 이런형식으로 되어 있기 때문에 [:-2] 를 하여 데이터 비교
            out.write(a)
        else:
            continue
    df.close()
    out.close()

main()