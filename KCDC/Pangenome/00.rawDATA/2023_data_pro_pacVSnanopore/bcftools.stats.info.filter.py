### after stats


import os,sys,glob


wDir = "/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/01.vcf.stat/"
file_lists = glob.glob(wDir+"*stats")

out = open("%sfilter.info.txt"%wDir,"w")
out.write("file\tsample\trecords\tno_ALTs\tSNPs\tMNPs\tindels\tothers\tmultiallelic_sites\tmultiallelic_SNP_sites\tts\ttv\tts_tv\tts_tv_1alt\n")

for i in file_lists:
    os.system("grep -v \"^#\" %s | head -11 > temp.filter.txt"%(i))
    tmp = open("temp.filter.txt","r")
    tmp1 = [a.replace("\n","") for a in tmp]
    #for row_index,j in enumerate(tmp1):
    for row_index,j in enumerate(tmp1):
        t2 = j.split("\t")
        if row_index == 0:
            out.write("%s"%j.split("\t")[2])
        elif t2[0] == "TSTV":
            out.write("\t%s\t%s\t%s\t%s"%(t2[2],t2[3],t2[4],t2[7]))
        else:
            out.write("\t%s"%(j.split("\t")[3]))
    out.write("\n")
        
os.system("rm temp.filter.txt")
out.close()

