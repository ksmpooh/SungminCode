# 2021 HLA long seq pro.
#HLA.Shortread.Seq.NIH19KT0748.align.sorted.dedup.bam
#samtools view -H $BAM | sed "s/Solid5500XL/Solid/" | samtools reheader - $BAM
import os,glob


#samtools addreplacerg -r "SM:${base_name}" -o ${base_name}.bam $file;
#samtools addreplacerg -r "ID:smkim123" -o addreplacerg.bam 2nd_Cell.264.bc1001--bc1001.bam

wDir = "/DATA/smkim/HLA_seq/MHC_PacBio_Rawdata/01.unmapped_bam/"
#out = "/BDATA/smkim/pangenome/01.revio_kchip/00.beforeMapping"
#KBA_ID	shortread_ID	Shortread_filename_R1	Shortread_filename_R2	longread_ID	Longread_filename	Longread_filePath
#NIH19KT0247	247	247_S99_L002_R1_001	247_S99_L002_R2_001	2020HLAseq001	1st_Cell.bc1001--bc1001	./MHC_PacBio_Rawdata/1st_Cell_CCS/1st_Cell.bc1001--bc1001


def main():
    dfs = glob.glob(wDir + "*.bam")
    #print(dfs)
    #2nd_Cell.734.bc1007--bc1007.bam
    refs = open("/DATA/smkim/HLA_seq/HLAseq.ID.table.txt","r")
    refs = [s.replace("\n","") for s in refs]
    ref_dic = {}
    for i in refs[1:]:
        line = i.split("\t")
        ref_dic.setdefault(line[1],line[0])
    
    for i in dfs:
        
        ##2nd_Cell.734.bc1007--bc1007.bam
        line = i.replace(wDir,"").split(".")
        idx = line[1]
        print("%s :%s"%(i,ref_dic[idx]))
        os.system("samtools view -H %s > tmp.header.txt"%i)
        old_header = open("./tmp.header.txt","r")
        old_header = [s for s in old_header]
        new_header = open("./tmp.new.header.txt","w")
        for j in old_header:
            if "@RG" in j:
                tmps = j.split("\t")
                for tmp in tmps:
                    if "SM:" in tmp:
                        old_SMtag = tmp
                        break
                new_header.write(j.replace(old_SMtag,"SM:%s"%(ref_dic[idx])))
            else:
                new_header.write(j)
        new_header.close()
        out = wDir + "HLA.Longread.Seq.%s.bam"%ref_dic[idx]
        os.system("samtools reheader tmp.new.header.txt %s > %s"%(i,out))


main()
        

    
#['@HD\tVN:1.5\tSO:unknown\tpb:3.0.1\n', '@RG\tID:6bbea546/0--0\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=101-789-500;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000;BarcodeFile=/denovo/workspace.bsy/bin/bin/sequel_report/Sequel_96_barcodes_v1.fasta;BarcodeHash=c7bb8c22c94e0982d3f06eb7665823d5;BarcodeCount=96;BarcodeMode=Symmetric;BarcodeQuality=Score\tLB:HLA_Pac_set3\tPU:m64224e_211223_091137\tSM:test123\tPM:SEQUELII\tBC:CACATATCAGAGTGCG\tCM:S/P4-C2/5.0-8M\n', '@PG\tID:ccs-6.2.0\tPN:ccs\tVN:6.2.0\tDS:Generate circular consensus sequences (ccs) from subreads.\tCL:ccs ccs -j 32 m64224e_211223_091137.subreads.bam m64224e_211223_091137.ccs.bam\n', '@PG\tID:lima\tVN:1.11.0 (commit v1.11.0)\tCL:lima m64224e_211223_091137.ccs.bam /denovo/workspace.bsy/bin/bin/sequel_report/Sequel_96_barcodes_v1.fasta Samples.bam --split-bam-named --peek-guess --same --dump-removed --min-score 0 --num-threads 32\n']

    


