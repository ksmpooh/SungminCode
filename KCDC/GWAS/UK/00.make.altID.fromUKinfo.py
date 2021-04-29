#/BDATA/smkim/UK/marker/MAFINFO
#ukb_mfi_chr13_v3.txt
#makeID
#1:10177_A_AC	rs367896724	10177	A	AC	0.40079	AC	0.467935
#1:10235_T_TA	rs540431307	10235	T	TA	0.000367353	TA	0.214688
#1:10352_T_TA	rs201106462	10352	T	TA	0.394625	TA	0.447895

import os

def main():
    outDir  = "/BDATA/smkim/UK/marker/MAFINFO/makeID/"
    for i in range(1,22+1):
        chr = str(i)
        os.system("awk '{print \"%s:\"$3\"_\"$4\"_\"$5\"\t\"$2}' ukb_mfi_chr%s_v3.txt > %sukb_mfi_chr%s_v3.newID.txt"%(chr,chr,outDir,chr))

main()
