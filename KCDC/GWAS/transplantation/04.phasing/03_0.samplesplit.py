####
# python 03_0.samplessplit.py INPUT.fam OUTPUT
####
import os, glob,sys



def main():
	fam = sys.argv[1]
	outDir = sys.argv[2]
	outDir = outDir + "/5KsplitSample/"
	os.system("mkdir %s"%outDir)

	df = open(fam,'r')
	size = len([a for a in open(fam,'r')])
	print("size : %i"%size)
	windowSize = 5000
	i = 1
	while 1:
		front = i
		tail = i+windowSize - 1
		if tail > size:
			tail = size
		print("front : %s, tail : %d"%(front,tail))
		out = "%s%s_%s.sampleID"%(outDir,front,tail)
		with open(out,'w') as output:
			for i in range(front,tail+1):
	                        output.write(df.readline().split(" ")[0] + "\n")
		if tail == size:
			break
		i = tail + 1
		print(i)

main()

