import os, glob
import multiprocessing as mp
from multiprocessing import Pool
import Queue

def doWork(x):
	print '%s' % (mp.current_process().name)
	os.system("sh " + x)
	return x

# read data
def readData(path, separator=' ', header=True):
	l = []
	with open(path, 'r') as f:
		if header is True:
			next(f)
		for line in f:
			l.append(line.strip().split(separator))
	return l

def main():
	#DATASET = readData("1000GP_Phase3_chr1.legend_10", ' ', True)

	folder = os.getcwd()
	DATASET = glob.glob('*.sh')
	
	#PROCESSES = max([1, mp.cpu_count()-1])
	#PROCESSES = max([1, (mp.cpu_count()/2)-4])
	PROCESSES = max([1, 25])
	print 'cpu_count() = %d\n' % PROCESSES

	pool = Pool(processes=PROCESSES)

	#with open('results.txt', 'w') as f:
	#	for result in pool.map(doWork, DATASET):
	#		f.write('%s\n' % result[0])

	pool.map(doWork, DATASET)
	
	pool.close()
	pool.join()

if __name__ == "__main__":
	main()

