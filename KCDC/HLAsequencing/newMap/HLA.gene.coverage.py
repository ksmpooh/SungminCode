'''
python HLA.gene.coverage.py [build] [wDir]

# build : hg19, hg38, hg38_alt
'''

import os,glob,sys

wDir = "/BDATA/smkim/HLA_seq/"
shDir = "/BDATA/smkim/HLA_seq/SCRIPTs/coverage/"
os.system("mkdir %s"%shDir)
#chr6_28510020_33480577:1-1000000
#0      1           2                       3   4       5       6       7       8   9                           10                          11  12      13      14      15
#812	NR_132323.1	chr6:28510020-33480577	+	1281885	1287787	1287787	1287787	3	1281885,1282416,1287172,	1282285,1282686,1287787,	0	HLA-V	none	none	-1,-1,-1,

genes = ["HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1"]

def make_target_ref(fileIn,target_contig):
    tmp = open(fileIn,'r')
    out = []
    while 1:
        line = tmp.readline().replace("\n","").split()
        if not line: break
        if (line[2] not in target_contig) or (line[12] not in genes):
            continue
        else:
            #[gene,contig,start,end]
            #[HLA-A,chr6_alt,3,1000]
            out.append(["%s,%s,%s,%s"%(line[12],line[2],line[4],line[5])])
    return out



        

# contig
def main():
    
    # build : hg19, hg38, hg38_alt
    build = sys.argv[1]
    wDir = sys.argv[2]
    if build == "hg19":
        ref = "/BDATA/smkim/HLA_seq/REF/IGV/hg19_ncbiRefSeq.sorted_onlyForHLA_forIGV.txt"
        target_contig= ["6_28477797_33448354"]
    elif build == "hg38":
        ref = "/BDATA/smkim/HLA_seq/REF/IGV/hg38_ncbiRefSeq_chr6_HLAregion_forIGV.txt"
        target_contig = ["chr6_28510020_33480577"]
    else:
        ref = "/BDATA/smkim/HLA_seq/REF/IGV/hg38_ncbiRefSeq_chr6_HLAregion_forIGV.txt"
        target_contig = ["chr6_28510020_33480577","chr6_GL000250v2_alt", "chr6_GL000251v2_alt", "chr6_GL000252v2_alt", "chr6_GL000253v2_alt", "chr6_GL000254v2_alt", "chr6_GL000255v2_alt", "chr6_GL000256v2_alt", "chr6_KI270758v1_alt"]
    print("ref : %s"%ref)
    print("build : %s"%(build))
    print("Working Directory : %s"%wDir)

    ref_df = make_target_ref(ref,target_contig)
    bams = glob.glob("%s/*.bam"%wDir)
    os.system("mkdir %s/coverage"%wDir)
    outDir = wDir + "/coverage"
    #shDir = 
    for bam in bams:
        for line in ref_df:
            line = line[0]
            gene, contig, start, end = line.split(",")
            tmp = "coverage_%s_%s"%(contig,gene)
            out = bam.replace(wDir,outDir).replace(".bam",".bam.%s"%tmp)
            shout = out.replace(outDir,shDir) + ".sh"
            #print(shout)
            shOut = open(shout,'w')
            shOut.write("samtools coverage -r %s:%s-%s %s > %s"%(contig,start,end,bam,out))
            shOut.close()

main()


    
    


