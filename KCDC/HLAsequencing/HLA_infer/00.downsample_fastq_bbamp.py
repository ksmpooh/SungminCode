### fastq downsampling
## python downsampling.py R1_input.fastq depth
## python downsampling.py R1.fastq 10
import os,glob,sys

'''
/BDATA/smkim/TOOLs/bbmap/reformat.sh \
in=/BDATA/smkim/HLA_seq/00.rawDATA/shortread/1003_S71_L002_R1_001_paired.fastq.gz \
out1=/BDATA/smkim/HLA_seq/downsampling/test/1003_S71_L002_R1_001_paired_down.10x.fastq.gz \
in2=/BDATA/smkim/HLA_seq/00.rawDATA/shortread/1003_S71_L002_R2_001_paired.fastq.gz \
out2=/BDATA/smkim/HLA_seq/downsampling/test/1003_S71_L002_R2_001_paired_down.10x.fastq.gz \
samplereadstarget=331370  \
overwrite=true
'''
#-R '@RG\tID:HWI\tSM:1003\tPL:ILLUMINA_down_10x\tLB:Novaseq6000' \
#bwa-mem2 mem -t 32 /BDATA/smkim/HLA_seq/REF/HLA.target.v2.fasta /BDATA/smkim/HLA_seq/downsampling/test/1003_S71_L002_R1_001_paired_down.10x.fastq.gz /BDATA/smkim/HLA_seq/downsampling/test/1003_S71_L002_R2_001_paired_down.10x.fastq.gz | samtools sort -o 1003.10x.test.bam


wDir = "/BDATA/smkim/HLA_seq/"
refDir = "/BDATA/smkim/HLA_seq/REF/"
short_dir = wDir + "00.rawDATA/shortread/"

downDir = wDir + "downsampling/"
shDir = downDir  + "SCRIPTs/"
os.system("mkdir %s"%shDir)
'''HLAseq.ID.table_v2.txt
KBA_ID	shortread_ID	Shortread_filename_R1	Shortread_filename_R2	longread_ID	Longread_filename	Longread_filePath	CELL	ID_check	OLD_ID	NEW_ID
NIH19KT0247	247	247_S99_L002_R1_001	247_S99_L002_R2_001	2020HLAseq001	KDCDP.2020HLAseq001.bc1001--bc1001.ccs	NA	1st_Cell_CCS	2020HLAseq001	NIH19KT0247	NIH19KT0247
'''


'''
def shortread_mapping_sh(ID_table,col_index):
    print("shortread mapping!")
    tool = "/BDATA/smkim/TOOLs/bwa-mem2/bwa-mem2"
    shDir2 = shDir + "01.shortread_mapping_%s/"%(theme)
    os.system("mkdir %s"%shDir2)
    for i in ID_table:
        i = i.split()
        #print(i)
        R1 = shortread_raw + i[2] + "_paired.fastq.gz"
        R2 = R1.replace("_R1_","_R2_")
        #out = shortread_mapping + "HLA.Shortread.Seq.%s.align.sorted.bam"%i[0]
        out = "%sHLA.Shortread.Seq.%s.trimmed.%s_align.sorted.bam"%(shortread_mapping,i[0],theme)
        with open(shDir2 + "%s_%s.sh"%(i[0],i[2]),"w") as shout:
            shout.write("%s mem -t 16 -R \"@RG\\tID:HWI\\tSM:%s\\tPL:ILLUMINA\\tLB:Novaseq7000\" %s %s %s | samtools sort -o %s\n"%(tool,i[0],ref,R1,R2,out))
            shout.write("samtools index %s\n"%(out))
'''
'''
#def shortread_mapping_sh(ID_table,col_index,):
def shortread_mapping_sh(ID_table,depth):
    print("shortread mapping!")
    tool = "/BDATA/smkim/TOOLs/bwa-mem2/bwa-mem2"
    #shDir2 = shDir + "01.shortread_mapping_%s/"%(theme)
    shDir = shDir + "02.shortread_mapping/"
    outDir = donwDir + "%s/01.mapping/"
    os.system("mkdir %s"%shDir)
    for i in ID_table:
        i = i.split()
        #print(i)
        R1 = shortread_raw + i[2] + "_paired.fastq.gz"
        R2 = R1.replace("_R1_","_R2_")
        #out = shortread_mapping + "HLA.Shortread.Seq.%s.align.sorted.bam"%i[0]
        out = "%sHLA.Shortread.Seq.%s.trimmed.%s_align.sorted.bam"%(shortread_mapping,i[0],theme)
        with open(shDir2 + "%s_%s.sh"%(i[0],i[2]),"w") as shout:
            shout.write("%s mem -t 16 -R \"@RG\\tID:HWI\\tSM:%s\\tPL:ILLUMINA\\tLB:Novaseq7000\" %s %s %s | samtools sort -o %s\n"%(tool,i[0],ref,R1,R2,out))
            shout.write("samtools index %s\n"%(out))
'''    


def downsampling(short_lists,depth):
    shDir = downDir  + "SCRIPTs/01.downsampling/"
    os.system("mkdir %s"%shDir)
    outDir = downDir + "%sX/"%depth
    os.system("mkdir %s"%(outDir))
    Tool = "/BDATA/smkim/TOOLs/bbmap/reformat.sh"
    region = 33448354-28477797
    basic_depth = int(region/101)
    target_depth = str(basic_depth * int(depth))
    for i in short_lists:
        in1 = i
        in2 = i.replace("_R1_","_R2_")
        out1 = in1.replace(short_dir,outDir).replace(".fastq.gz",".%sX.fastq.gz"%depth)
        out2 = in2.replace(short_dir,outDir).replace(".fastq.gz",".%sX.fastq.gz"%depth)
        with open(shDir + "%sX_%s.sh"%(depth,i.replace(short_dir,"")),"w") as shout:
            shout.write("%s in=%s out1=%s in2=%s out2=%s samplereadstarget=%s overwrite=true"%(Tool,in1,out1,in2,out2,target_depth))
            
        



def main():
    ref = "/BDATA/smkim/HLA_seq/REF/HLA.target.v2.fasta"
    refs = open(wDir + "HLAseq.ID.table_v2.txt","r")
    refs = [s.replace("\n","") for s in refs]
    header = refs.pop(0)
    short_R1_list = glob.glob(short_dir + "*_R1_*.fastq.gz")

    for depth in ["10","20","30","50","100","200"]:
        downsampling(short_R1_list,depth)
        # mapping after trimming
        


main()





