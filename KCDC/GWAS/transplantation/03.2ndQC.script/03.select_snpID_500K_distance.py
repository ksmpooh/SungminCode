import os,random

wdir = "/ADATA/smkim/JG/03.QC_2nd/ethicalPCA/"

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
bim = fileRead(wdir + "merge.bim")



windowSize = 500000

before_chrom = '0'
#matrix = []
#profile = []
index = []

for i,line in enumerate(bim):
	chrom,ID,trash,position,ref,alt = line.split('\t')
	if before_chrom != chrom:
		front = int(position)
		index.append(i)
		before_chrom = chrom
		continue
	if int(position) >= front + windowSize:
		#profile = [chrom,front,int(position),front_index,i-1]
		#index.append(random.randrange(front_index,i))
		index.append(i)
		front_index = i
		front = int(position)
		#matrix.append(profile)


for i in index:
	randomID_list.write(bim[i].split("\t")[1] + "\n")
randomID_list.close()

#fileWrite(matrix,wdir+"chr_front_last_frontIndex_last_index.txt")

