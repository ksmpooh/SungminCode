import os,glob

wDir = "/BDATA/smkim/JG_2020/08.asso/OUTPUTs/"

mergeDir = wDir + "01_1.merge/"
plotDir = mergeDir + "plot/"
os.system("mkdir " + plotDir)
#plotsh = "/BDATA/smkim/JG/08.asso/SCRIPTs/03.ploting.R"
#Three args : [input.txt] [output]

def main():
    #phenos = ["HT","HT_T2D","sub_Total","T2D","Total"]
    dfs = glob.glob(mergeDir + "ESRD*.txt")
    for df in dfs:
        print("df : " + df)
        datain = df
        os.system("grep -v ^# %s |awk '{print $4,$1,$2,$9}' > %s"%(datain, datain.replace(".txt",".manhattan.form.txt"))) 
#        plotout = datain.replace(mergeDir,plotDir).replace(".txt","")
#        datain = datain.replace(".txt",".manhattan.form.txt")
#        os.system("Rscript --vanilla %s %s %s"%(plotsh,datain,plotout))

main()








