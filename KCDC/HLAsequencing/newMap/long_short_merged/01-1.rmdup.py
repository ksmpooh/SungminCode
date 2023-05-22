import os,glob,sys

## python 01.rmdup.py [bamdir] [outDir]
bamDir = sys.argv[1]
outDir = sys.argv[2]
os.system("mkdir %s"%outDir)
refDir = "/DATA/smkim/HLA_seq/REF/"

## input sorted bam -> dup.bam

def mark_duplicated():
    print("Mark_duplicated")
    dfs = glob.glob(bamDir + "*.sort.bam")
    shDir2 = "./01-1.shortread_markDup/"
    os.system("mkdir %s"%shDir2)
    for i in dfs:
        i = i.replace(bamDir,"")
        out = i.replace(".bam",".dedup.bam")
        with open(shDir2 + out.replace(".bam",".sh"),"w") as shout:
            shout.write("sudo docker run -v \"%s\":\"/input\" -v \"%s\":\"/ref\" broadinstitute/picard java -jar /usr/picard/picard.jar MarkDuplicates I=/input/%s O=/input/%s METRICS_FILE=duplicates REMOVE_DUPLICATES=True CREATE_INDEX=True\n"%(bamDir,refDir,i,out))
            out = outDir + "/" +  out
            shout.write("samtools index %s\n"%(out))
            #shout.write("samtools stats %s > %s.stats\n"%(out,out))


def main():
#    refs = open(wDir + "HLAseq.ID.table.txt","r")
#    refs = [s.replace("\n","") for s in refs]

#    header = refs.pop(0)
    ##shortread_mapping_sh(refs,2)
    mark_duplicated()
    #bqsr()
    #deepvariant_calling()
    #gatk_calling()

main()