####
# python 03_0.samplessplit.py INPUT.fam OUTPUT
####
import os, glob,sys



def main():
	fam = sys.argv[1]
	outDir = sys.argv[2]
	#size = sys.argv[3]
	outDir = outDir + "/5KsplitSample/"
	os.system("mkdir %s"%outDir)
	
	df = open(fam,'r')
	size = len([a for a in open(fam,'r')])
	print("size : %i"%size)
	windowSize = 5000
	i = 1
	while 1:
#	for i in range(1,len(fam)):
		front = i
		tail = i+windowSize - 1
		print("front : %d, tail : %d"%(front,tail))
		if tail > size:
			tail = size
		out = "%s%s_%s.sampleID"%(outDir,front,tail)
		with open(out,'w') as output:
			for i in range(front,tail+1):
	                        output.write(df.readline().split(" ")[0] + "\n")
		if tail == size:
			break
		i = tail + 1
		print(i)

main()

