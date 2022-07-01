import os, glob


wDir = "/LaCie2/KOTRY/99.open/forOpen/"
outDir = "/BDATA/smkim/JG/99.openforKOTRY/02.Imputed.VCF/"
shDir = outDir + "SCRIPTs/"
os.system("mkdir %s"%shDir)

organs = ["KD","KR","LR","LD"]
QCs = ["discovery","replication"]

#/LaCie2/KOTRY/99.open/forOpen/KD/discovery/02.Imputed.VCF

#bcftools concat -Oz >
def main():
    for organ in organs:
        print(organ)
        for QC in QCs:
            out = open(shDir + "%s_%s.merge.sh"%(organ,QC),"w") 
            out.write("bcftools concat -Oz ")
            for chrom in range(1,22+1):
                print(str(chrom))
                tmp = glob.glob("/LaCie2/KOTRY/99.open/forOpen/%s/%s/02.Imputed.VCF/*chr%s_*gz"%(organ,QC,str(chrom)))
                #print(tmp)
                out.write("%s "%(tmp.pop()))
                #KBA.KOTRY.LR.replication.QCed.PLINK
            out.write("-Oz > %sKBA.KOTRY.%s.%s.IMPUTE4_IMPUTED.filterMAF0.01_INFO0.08.vcf.gz\n"%(outDir,organ,QC))


#main()
import itertools
def make_new_header(ref_dic,VCFin):
    os.system("bcftools view -h %s > tmp.header.txt"%(VCFin))
    header_in = open("tmp.header.txt","r")
    headers = [s for s in header_in]
    header_out = open("new.header.txt","w")
    count = 0 
    for header in headers:
        count = count + 1
        if count == len(headers):
            i = header.split()
            for j in i[0:7+1]:
                header_out.write("%s\t"%j)
            header_out.write("FORMAT")
            for j in i[9:]:
                header_out.write("\t%s"%ref_dic[j])
            header_out.write("\n")
            return 
        header_out.write("%s"%header)


def header_change():
    ref = open("/BDATA/smkim/JG/00.rawData/JG.update.ID.NIHtobCODE.txt","r")
#KBA_ID KBA_ID bCODE bCODE
#NIH19KT0001 NIH19KT0001 02903093DNA01102 02903093DNA01102
    ref_list = [s.replace("\n","") for s in ref]
    print(ref_list[1:5])
    ref_dic = {}
    for ref in ref_list[1:]:
        oldID,oldID2,newID,newID2 = ref.split()
        ref_dic[oldID] = newID
#    print(dict(itertools.islice(ref_dic.items(), 5)))
    make_new_header(ref_dic,"KBA.KOTRY.KD.discovery.IMPUTE4_IMPUTED.filterMAF0.01_INFO0.08.vcf.gz")



#header_change()

#bcftools view -h input.vcf >  header.txt


