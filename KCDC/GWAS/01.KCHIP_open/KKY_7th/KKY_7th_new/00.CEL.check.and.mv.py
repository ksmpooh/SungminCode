import glob,os

wdir = "/BDATA/smkim/KKY_7th/00.rawData/"

CEL = wdir + "CEL/2020_pro_7th_Teragen/"
out = CEL + "out/"

def main():
    refin = open(wdir + "2020.all.teragen.KKY.info_20220311_only6th.txt","r")
    refs = [s.replace("\n","").split()[0] for s in refin]
    print(refs[0:5])
    for ref in refs:
        mv = glob.glob(CEL + "*%s*.CEL"%(ref))
        os.system("mv %s %s"%(mv.pop(),out))


main()