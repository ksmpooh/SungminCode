### michigan vs original
### HLA imputation result preprocessing
# python 02.extract.HLA.py input.vcf.gz output.txt

## after KIS KHU
# zcat < chr6.info.gz | grep -v rs |grep -v chr |grep -v gg |grep -v SNPS |grep -v AA | awk '{print $1}' > hla.type.list.txt
#  cat < chr6.dose.vcf.gz | grep -v rs |grep -v chr6_ |grep -v gg |grep -v SNPS |grep -v AA | bgzip -c > chr6.dose.only.HLA.vcf.gz

'''
HLA_A*01    NIH19KT0019_NIH19KT0019=0	NIH19KT0013_NIH19KT0013=0	NIH19KT0008_NIH19KT0008=0
HLA_A*01:01:01:01   NIH19KT0019_NIH19KT0019=0	NIH19KT0013_NIH19KT0013=0	NIH19KT0008_NIH19KT0008=0
HLA_A*01:01:13	NIH19KT0019_NIH19KT0019=0	NIH19KT0013_NIH19KT0013=0	NIH19KT0008_NIH19KT0008=0
HLA_A*01:02 NIH19KT0019_NIH19KT0019=0	NIH19KT0013_NIH19KT0013=0	NIH19KT0008_NIH19KT0008=0
HLA_A*01:12 NIH19KT0019_NIH19KT0019=0	NIH19KT0013_NIH19KT0013=0	NIH19KT0008_NIH19KT0008=0
HLA_A*01:136    NIH19KT0019_NIH19KT0019=0	NIH19KT0013_NIH19KT0013=0	NIH19KT0008_NIH19KT0008=0
HLA_A*02    NIH19KT0019_NIH19KT0019=0.001	NIH19KT0013_NIH19KT0013=0.004	NIH19KT0008_NIH19KT0008=0
HLA_A*02:01:01:01   NIH19KT0019_NIH19KT0019=0.001	NIH19KT0013_NIH19KT0013=0.003	NIH19KT0008_NIH19KT0008=0
HLA_A*02:01:02  NIH19KT0019_NIH19KT0019=0	NIH19KT0013_NIH19KT0013=0	NIH19KT0008_NIH19KT0008=0
HLA_A*02:01:04  NIH19KT0019_NIH19KT0019=0	NIH19KT0013_NIH19KT0013=0	NIH19KT0008_NIH19KT0008=0

to

ID      HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2
NIH19KT0019     A*33    A*24    B*07    B*58
NIH19KT0013     A*24    A*02    B*40    B*15
NIH19KT0008     A*03    A*26    B*44    B*15
NIH19KT0254     A*02    A*33    B*07    B*46
NIH19KT0264     A*03    A*02    B*44    B*37
NIH19KT0410     A*02    A*68    B*27    B*38
NIH19KT0415     A*02    A*01    B*15    B*37
NIH19KT0747     A*26    A*24    B*15    B*40

'''

## python vcfin outfile.txt #digit(2/4)
import os,glob,sys

vcfIn = sys.argv[1]
#vcfIn = "/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/michigan/chr6.dose.vcf.gz"
out_path = sys.argv[2]
#out_path = "/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/michigan/python_test.txt"
#digit = sys.argv[3]

if "1KGP" in vcfIn:
    target_genes = ["HLA_A","HLA_B","HLA_C","HLA_DQB1","HLA_DRB1"]
else:
    target_genes = ["HLA_A","HLA_B","HLA_C","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1","HLA_DRB1"]


#target_genes = ["HLA_A","HLA_B","HLA_C","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1","HLA_DRB1"]
#target_genes = ["HLA_A","HLA_B","HLA_C","HLA_DQB1","HLA_DRB1"]



'''

new = [s for s in header.split("\t")]
c = len(input_list[0].split("\t"))
print(input_list[0].split("\t"))
print(c)
for i in input_list:
    t = i.split()
    for j in range(0,c+1):
        if j ==0:
            new[j] = new[j] +"\t" + t[j]#.split("=")
        else:
            new[j] = new[j] +"\t" + t[j].split("=")[1]


for i in new:
    print(i)
'''

def t_dataframe(df,gene):
    header = "ID\t"
    #print(df[0])
    tmp = df[0].strip().split()
    # ID_ID=0 나누기 -> ID header
    #tmp = [s.split("_")[0] for s in tmp[1:]]
    tmp = [s for s in tmp[1:]]
    '''
    if "_" in tmp[1:]:
        tmp = [s.split("_")[0] for s in tmp[1:]]
    else:
        tmp = [s.split("=")[0] for s in tmp[1:]]
    '''
    header = header + '\t'.join(tmp)
    t_df = [s for s in header.split()]
    #print(t_df)
    c = len(df[0].split())
    for i in df:
        t = i.split()
        for j in range(0,c):
            if j == 0:
                t_df[j] = t_df[j] + "\t" + t[j]
            else:
                t_df[j] = t_df[j] + "\t" + t[j].split("=")[1]
    return t_df


 #   out = t_df[0]
#    out = [s.split() for s in t_df[1:]]
def digit_split(df,gene): #2 or other
    print("%s"%(gene))
    td =[]
    fd =[]
    #header_o = df[0]
    #header = "ID"

    for i in df:
        tmp = i.split()
        g1,allele = tmp[0].split("*")
        tmp = [tmp[0]]
        tmp = [s for s in tmp]
        if len(allele.split(":")) == 1:
            td.append(i)
        else:
            fd.append(i)
    return t_dataframe(td,gene),t_dataframe(fd,gene)


def HLAtyping_select(df,gene):
    out = []
    header = df[0].split()
    
    #out.append("ID\t%s_1\t%s_2"%(gene,gene))
    out.append("ID\t%s.1\t%s.2"%(gene,gene))

    for i in df[1:]:
        tmp = i.split()
        sample_ID = tmp[0]
        tmp[0] = "-1"
        type_list = [float(s) for s in tmp]
        max1_idx = type_list.index(max(type_list))
        type_list[max1_idx] = -1
        max2_idx = type_list.index(max(type_list))
        t1 = header[max1_idx].replace("HLA_","")
        t2 = header[max2_idx].replace("HLA_","")
        out.append("%s\t%s\t%s"%(sample_ID,t1,t2))

    return out
    

def df_merge(out,df,gene):
    if gene == "HLA_A":
        return df
    for i in range(0,len(df)):
        out[i] = out[i] + "\t" + "\t".join(df[i].split("\t")[1:])
    return out
        





def main():
    print("vcf to query")
    #bcftools query -f '%ID\t%REF\t%ALT[\t%SAMPLE=%DS]\n' chr6.dose.vcf.gz |grep HLA > test.txt
    td_a = []
    fd_a = []
    for gene in target_genes:
    #for gene in target_genes[0:2]:
        print("Target gene : %s"%gene)
        tmp_path = vcfIn + ".tmp"
        os.system("bcftools query -f '%%ID[\t%%SAMPLE=%%DS]\n' %s | grep %s > %s"%(vcfIn,gene,tmp_path))
        #os.system("bcftools query -f \'%%ID[\t%%SAMPLE=%%DS]\n\' chr6.dose.vcf.gz | grep HLA_A > %s_tmp"%("hi"))
        #tmp = open("tmp","r")
        tmp = open(tmp_path,"r")
        tmp = tmp.readlines()
        print("DATA preprocessing...: %s"%gene)
        td,fd  = digit_split(tmp,gene)
        td = HLAtyping_select(td,gene)
        fd = HLAtyping_select(fd,gene)
        #print(td)
        #print(t_dataframe(td,gene))
        print("DATA merge...%s"%gene)
        td_a = df_merge(td_a,td,gene)
        fd_a = df_merge(fd_a,fd,gene)
    #os.system("rm %s"%tmp_path)

    td_out = open(out_path.replace(".txt","_td.txt"),'w')
    fd_out = open(out_path.replace(".txt","_fd.txt"),'w')
    for t in td_a:
        td_out.write("%s\n"%t)
        #td_out.write("%s\n"%"\t".join(t))
        #fd_out.write("%s\n"%"\t".join(f))
    for f in fd_a:
        fd_out.write("%s\n"%f)

    td_out.close()
    fd_out.close()
    


main()
