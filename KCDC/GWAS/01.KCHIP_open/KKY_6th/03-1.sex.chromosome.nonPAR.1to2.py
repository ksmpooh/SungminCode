#cd /BDATA/smkim/KKY_6th/03.2ndQC/sex

import os


dfs = open("test.raw","r")
out = open("test.grep1.txt","w")
dfs = [s.replace("\n","") for s in dfs]
out.write(dfs.pop(0))
out.write("\n")
for df in dfs:
     lines = df.split()
     line = lines[6:]
     if lines[4] == "1" and line.count("1") != 0:
          out.write(df)
          out.write("\n")
          
out.close()






dfs = open("test.raw","r")
#out = open("test.grep1.txt","w")
out = open("test_1to2.txt","w")
dfs = [s.replace("\n","") for s in dfs]
out.write(dfs.pop(0))
out.write("\n")
for df in dfs:
    lines = df.split()
    line = lines[6:]
    if lines[4] == "1" and line.count("1") != 0:
        idx = ' '.join(lines[0:6])
        line = ' '.join(line).replace("1","2")
        out.write(idx)
        out.write(" ")
        out.write(line)
        out.write("\n")
    else:
        out.write(df)
        out.write("\n")


out.close()




import os


dfs = open("KKY.6th.preQC.onlySEX_nonPAR_fil_raw.raw","r")
out = open("KKY.6th.preQC.onlySEX_nonPAR_fil_raw_grep1.txt","w")
dfs = [s.replace("\n","") for s in dfs]
out.write(dfs.pop(0))
out.write("\n")
for df in dfs:
     lines = df.split()
     line = lines[6:]
     if lines[4] == "1" and line.count("1") != 0:
          out.write(df)
          out.write("\n")

out.close()









import os

inDATA = "KKY.6th.preQC.onlySEX_nonPAR_recode.ped"
dfs = open(inDATA,"r")
#out = open("test.grep1.txt","w")

out = open(inDATA.replace(".ped","_hettohomo.ped"),"w")
dfs = [s.replace("\n","") for s in dfs]
out.write(dfs.pop(0))
out.write("\n")
for df in dfs:
    lines = df.split()
    line = lines[6:]
    if lines[4] == "1":
        out.write(' '.join(lines[0:6]) + " ")
        while 1:
            if len(line) == 0:
                break
            a1 = line.pop(0)
            a2 = line.pop(0)
            if a1 != a2:
                a1 = "0"
                a2 = "0"
            out.write(a1 +" "+a2+" ")
        out.write("\n")
    else:
        out.write(df)
        out.write("\n")


out.close()