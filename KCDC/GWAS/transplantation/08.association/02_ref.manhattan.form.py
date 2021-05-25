#python 01.manhattan.form.py [datain : all] [trait]
#python 01.manhattan.form.py /BDATA/smkim/JG/08.asso/OUTPUTs/01_1.merge/JG.KR.NODAT.epacts.b.firth.result.txt JG.KR.NODAT.manhattan

import os, sys,glob


datain = sys.argv[1]
plotout = sys.argv[2]
#awk '(0.99>= $7 >=0.01 || 0.99 >= $6 >= 0.01) && ($15>=0.001){print $6,$7,$15}'
#awk '(0.99>= $7 >=0.01 || 0.99 >= $6 >= 0.01) && ($15>=0.001)
#{split($1,a,":");split(a[2],b,"_"); print a[1]"\t"b[1]"\t"$10}' AG_META_1.txt > ../plot/AG_META_KBA.forplot.txt
#inDir = "/DATA/smkim/Anthro/RESULTs/KBAmeta/"
#inDir = "/DATA/smkim/Anthro/RESULTs/BBJmeta/"

#outDir = "/DATA/smkim/Anthro/RESULTs/plot/"


#traits = glob.glob(inDir + "*.txt")

#for trait in traits:
#       os.system(awk '{split($1,a,":"); split(a[2],b,"_"); print a[1]"\t"b[1]"\t"$10}' AG_META_1.txt |

os.system("grep -v \"#\" %s | awk '{print $4,$1,$2,$9}' > %s"%(datain, datain.replace(".txt",".manhattan.form.txt")))

#Rscript --vanilla 03.ploting.R WBC_META_BBJ.with.KBA.forplot.txt WBC_KBAwithBBJ_meta
datain = datain.replace(".txt",".manhattan.form.txt")
os.system("Rscript --vanilla 03.ploting.R %s %s"%(datain,plotout))








