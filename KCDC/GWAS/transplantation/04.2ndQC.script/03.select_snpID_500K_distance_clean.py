
import os,random

wdir = "~/2ndQC/"

def fileRead(fileIn):
	f = open(fileIn,'r')
	inData = [r.replace("\n","") for r in f]
	return inData
def fileWrite(data,out):
	outFile = open(out,'w')
	for i in data:
		outFile.write(str(i)+ '\n')
	outFile.close()

randomID_list = open(wdir + "500K_distance_SNP_list.txt","w")
bim = fileRead(wdir + "ethnic.merge.bim")



windowSize = 500000

before_chrom = '0'
index = []

for i,line in enumerate(bim):
	chrom,ID,trash,position,ref,alt = line.split('\t')
	if before_chrom != chrom:
		front = int(position)
		index.append(i)
		before_chrom = chrom
		continue
	if int(position) >= front + windowSize:
		index.append(i)
		front_index = i
		front = int(position)


for i in index:
	randomID_list.write(bim[i].split("\t")[1] + "\n")
randomID_list.close()




