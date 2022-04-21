import os,glob,sys


cel = "/backup/2020_pro/"
#out_cel = "/backup/2020_pro/preg/"
#out_cel = "/backup/2020_pro/CAD/"
out_cel = "/backup/2020_pro/KKY/"

#refIn data
#NIH20G2818314	preg
#NIH20G2399419	preg
#NIH20G2733498	preg

def main():
    #refIn = open("/backup/smkim/preg_list.txt","r")
    #refIn = open("/backup/smkim/CAD_list.txt","r")
    refIn = open("/backup/smkim/KKY_list.txt","r")
    refs = [s.replace("\n","").split()[0] for s in refIn]
    print("CEL : %s"%str(len(refs)))
#    CELs = glob.glob(cel + "*CEL")
    count = 0
    for ref in refs:
        CEL = glob.glob(cel + "*%s*"%(ref))
        os.system("mv %s %s"%(CEL.pop(),out_cel))
        count = count + 1
        if count%100 == 0:
            print("count : "+str(count))
    print("count : "+str(count))

main()


# preg : 3623
# CAD : 3182
# KKY : 4525