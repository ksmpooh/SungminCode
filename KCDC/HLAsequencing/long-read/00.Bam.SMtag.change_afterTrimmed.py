# 20221103 trimmging 후에 메핑완료된 파일에 대하여 아이디 처리 및 MQ > 20 QC 
# 추가로 long-read 아이디 이상한거까지 갗이 처리 2<>4 3<>5
import os,glob

wDir = "/BDATA/smkim/HLA_seq/longread/"

inDir = wDir + "01.unmapped/trimmed/"
outDir = wDir + "02.mapped/"

#KBA_ID	shortread_ID	Shortread_filename_R1	Shortread_filename_R2	longread_ID	Longread_filename	Longread_filePath	CELL	ID_check	OLD_ID	NEW_ID
#NIH19KT0247	247	247_S99_L002_R1_001	247_S99_L002_R2_001	2020HLAseq001	KDCDP.2020HLAseq001.bc1001--bc1001.ccs	NA	1st_Cell_CCS	2020HLAseq001	NIH19KT0247	NIH19KT0247
'''
0 KBA_ID
1 shortread_ID
2 Shortread_filename_R1
3 Shortread_filename_R2
4 longread_ID
5 Longread_filename
6 Longread_filePath
7 CELL
8 ID_check
9 OLD_ID
10 NEW_ID
'''
'''
samtools view -H 5th_Cell.3832.bc1008--bc1008_trimmed_toBAM_mapped.bam > header.txt
@RG/tID:myid/tSM:mysample # 해더에 추가
samtools reheader header.txt 5th_Cell.3832.bc1008--bc1008_trimmed_toBAM_mapped.bam | samtools view -h -q 20 -o header.change.test.bam
'''        

#2nd_Cell.724.bc1002--bc1002_trimmed_toBAM_mapped.bam
#4th_Cell.bc1010--bc1010_trimmed_toBAM_mapped.bam
#5th_Cell.3832.bc1008--bc1008_trimmed_toBAM_mapped.bam
#KCDCP.2020HLAseq002.bc1002--bc1002.ccs_trimmed_toBAM_mapped.bam
def main():
    dfs = glob.glob(inDir + "*_mapped.bam")
    refs = open("/BDATA/smkim/HLA_seq/HLAseq.ID.table_v2.txt")
    refs = [s.replace("\n","") for s in refs]
    ref_dic = {}
    for i in refs[1:]:
        line = i.split("\t")
        ref_dic.setdefault(line[5],line[10])
    #print(refs)
    for i in dfs:
        line = i.replace(inDir,"").split(".")
        check = line[0]
        if check not in ["KCDCP","4th_Cell"]:
            line.pop(1)
        idx = '.'.join(line).replace("_trimmed_toBAM_mapped.bam","")
        print(idx)
        os.system("samtools view -H %s > tmp.header.txt"%i)
        old_header = open("tmp.header.txt","r")
        new_header = open("new.header.txt","w")
        old_df = old_header.readlines()
        for j in old_df:
            new_header.write(j)
        #print(old_df)
        #new_header.write("\n".join(old_df))
        new_header.write("@RG\tID:Pacbio_HLA_Hifi_SequalII SM:%s\n"%(ref_dic[idx]))
        new_header.close()
        out = outDir + "HLA.Longread.Seq.%s.trimmed.align.Q20.bam"%(ref_dic[idx])
        if os.path.isfile(out):
            continue
        os.system("samtools reheader new.header.txt %s | samtools view -h -q 20 -o %s"%(i,out))

#samtools reheader header.txt 5th_Cell.3832.bc1008--bc1008_trimmed_toBAM_mapped.bam | samtools view -h -q 20 -o header.change.test.bam
#HLA.Longread.Seq.%s.bam"%ref_dic[idx]
#KCDCP.2020HLAseq002.bc1002--bc1002.ccs_trimmed_toBAM_mapped.bam
main()
