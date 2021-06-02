
wDir = "/Users/ksmpooh/Desktop/KCDC/HLAimputation/patent/marker_list/"

def main():
    df = open(wDir + "KBA.ref.allele.for.HLAregion.txt","r")
    out = open(wDir + "forSNPnexus.input.v2.txt","w")
    header = df.readline()
    while 1:
        a = df.readline()
        if not a:
            break
        pos,ref,chr,KBA_ID,x,A1,A2 = a.split()
        if ref == A1:
            out.write("Chromosome\t6\t%s\t%s\t%s\t1\n"%(pos,A1,A2))
        elif ref == A2:
            out.write("Chromosome\t6\t%s\t%s\t%s\t1\n"%(pos,A2,A1))
        else:
            print("check : "+a)
    out.close()

main()
        