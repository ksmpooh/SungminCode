'''
id position a0 a1 all.aaf
10:60723:T:C 60723 T C 0.000310097
10:61064:T:C 61064 T C 0.000186058
10:61141:G:A 61141 G A 0.00130241
10:61420:G:A 61420 G A 0.000248077
10:61420:G:C 61420 G C 0.000744232
10:61512:G:A 61512 G A 0.000806252
10:61653:C:T 61653 C T 0.00136443
10:61657:T:C 61657 T C 0.00707021
'''

import os, glob,sys

filein = "chr10_wgs8k_imputationPanel.legend.gz"
fileout = "test.txt"
out = open(fileout,'w')
os.system("zcat %s > tmp"%filein)

df = open("tmp","r")
lines = df.readlines()
header = lines.pop(0)

start=0
end = 0
pre_pos = 0
chunk=5000000
#for line in df:
while 1:
    if len(lines) == 0:
        out.write("%s %s %s\n"%(chrom,start,pos))
        break
    line = lines.pop(0)
    line = line.split()
    chrom,pos,ref,alt=line[0].split(":")
    
    if start == 0:
        start = pre_pos
    elif (int(start) + int(chunk) <= int(pos)) or (int(pos) - int(pre_pos) > 1000000):
        print("prepos: %s, pos: %s"%(pre_pos,pos))
        #end = pos
        out.write("%s %s %s\n"%(chrom,start,pre_pos))
        start = 0
    else:
        pass
    pre_pos = pos
    


