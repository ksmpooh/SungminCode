import os,glob



def main():
    outDir = "/BDATA/smkim/HLAtyping/01.HLAgeneMerge/"
    years = ["2019","2020"]
    for year in years:
        if year == '2019':
            refs = "/BDATA/smkim/HLAtyping/00.rawDATA/HLAtypingIDmatching/2019.HLAtypingIDfolder.txt"
            wDir = "/BDATA/smkim/HLAtyping/00.rawDATA/2019/CDC_dbSNP_files/"
        else:
            refs = "/BDATA/smkim/HLAtyping/00.rawDATA/HLAtypingIDmatching/2020.HLAtypingIDfolder.txt"
            #wDir = "/Users/ksmpooh/Desktop/KCDC/일반용역/장기이식/2020_NGgene결과/2.CDC_dbSNP_files/"2.2020_265samples_11locus_variants
            wDir = "/BDATA/smkim/HLAtyping/00.rawDATA/2020/2020_265samples_11locus_엔젠바이오/2.2020_265samples_11locus_variants/"
        refs = open(refs,"r")    
        refs = [s.replace("\n","") for s in refs]
        print(refs[1:5])
        ref_dic = {}
        for i in refs[1:]:
            line = i.split("\t")
            ref_dic.setdefault(line[2],line[0])
        
        dfs = glob.glob(wDir + "*")
        print(dfs[0:5])
        for df in dfs:
            if ref_dic.get(df.replace(wDir,"")) == None:
                continue
            pre_vcf_list = glob.glob(df + "/*.vcf")
            post_vcf_list = [s for s in pre_vcf_list if os.path.getsize(s)!=0]  #os.path.getsize()
            output = outDir + "HLAtyping.VCF.%s.vcf"%(ref_dic[df.replace(wDir,"")])
            #print("bcftools concat %s/*vcf > %s"%(df,output))
            #os.system("bcftools concat %s/*vcf > %s"%(df,output))
            os.system("bcftools concat %s > %s"%(' '.join(post_vcf_list),output))
            #print(df)
            #print(post_vcf_list)
            #print(pre_vcf_list)
        

#main()


def VCF_name_change():
    wDir = "/BDATA/smkim/HLAtyping/"
    inDir = wDir + "01.HLAgeneMerge/"
    outDir = wDir + "02.IDchange/"
    vcfs = glob.glob(inDir + "*vcf")
    for vcf in vcfs:
        #unknown
        #HLAtyping.VCF.NIH19KT0252.vcf
        newID = vcf.replace(inDir,"").replace("HLAtyping.VCF.","").replace(".vcf","")
        tmpOut = open(wDir + "tmp.txt","w")
        tmpOut.write("unknown\t%s"%newID)
        #os.system("bcftools reheader --samples %s %s | bcftools sort - -Oz -o %s.gz"%(wDir+"tmp.txt",vcf,vcf.replace(inDir,outDir)))
        os.system("bcftools reheader --samples %s %s | bcftools sort -Oz -o %s.gz"%(wDir+"tmp.txt",vcf,vcf.replace(inDir,outDir)))
        #os.system("bcftools sort %s | bcftools reheader --samples %s -Oz -o %s"%(vcf,wDir+"tmp.txt",vcf.replace(inDir,outDir)))

VCF_name_change()
