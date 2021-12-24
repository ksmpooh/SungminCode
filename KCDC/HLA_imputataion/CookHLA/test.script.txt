 python -m MakeGeneticMap \
        -i example/1958BC.hg19 \
        -hg 19 \
        -ref 1000G_REF/1000G_REF.EUR.chr6.hg18.29mb-34mb.inT1DGC \
        -o MyAGM/1958BC+1000G_REF.EUR


 python CookHLA.py \
    -i example/1958BC.hg19 \
    -hg 19 \
    -o MyHLAImputation/1958BC+HM_CEU_REF \
    -ref example/HM_CEU_REF \
    -gm example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB \
    -ae example/AGM.1958BC+HM_CEU_REF.aver.erate \
    -mem 2g \
    # -mp 2   # The number of available cores for Multiprocessing.
    # -bgl4   # If you want to use Beagle4 instead of Beagle5.



python CookHLA.py \
    -i example/1958BC.hg19 \
    -hg 19 \
    -o MyHLAImputation/1958BC+HM_CEU_REF \
    -ref 1000G_REF/1000G_REF.EUR.chr6.hg18.29mb-34mb.inT1DGC
    -gm MyAGM/1958BC+1000G_REF.EUR.mach_step.avg.clpsB
    -ae MyAGM/1958BC+1000G_REF.EUR.aver.erate
    -mem 4g \
    -mp 8




python -m MakeGeneticMap \
-i JG/hg18/JG.QCed.HLA_rmAmbiguous \
-hg 18 -ref han/HAN -o JG.map/JG.map



python CookHLA.py \
    -i JG/hg18/JG.QCed.HLA_rmAmbiguous \
    -hg 18 \
    -o ./JG.imp/JG.Han.HG18.HLAimp \
    -ref ./han/HAN \
    -gm ./JG.map.mach_step.avg.clpsB \
    -ae ./JG.map.aver.erate \
    -mem 10g \
    -mp 8


python CookHLA.py     -i JG/hg18/JG.QCed.HLA     -hg 18     -o ./JG.imp/JG.Han.HG18.HLAimp     -ref ./han/HAN     -gm ./JG.map/JG.map.mach_step.avg.clpsB     -ae ./JG.map/JG.map.aver.erate     -mem 10g     -mp 8


python -m MakeGeneticMap \
-i /DATA/smkim/cookHLA/ESRD/JG.KR.merge.chr6 \
-hg 19 -ref han/HAN -o /DATA/smkim/cookHLA/ESRD/map/ESRD.map


python CookHLA.py \
    -i JG/hg18/JG.QCed.HLA_rmAmbiguous \
    -hg 18 \
    -o ./JG.imp/JG.Han.HG18.HLAimp \
    -ref ./han/HAN \
    -gm ./JG.map.mach_step.avg.clpsB \
    -ae ./JG.map.aver.erate \
    -mem 10g \
    -mp 8
