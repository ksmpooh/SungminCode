'''
Created on 2017. 6. 1.

@author: HMY
'''

import os, sys, glob, time

def shRUN(logDir):
        folder = os.getcwd()
        fileList = glob.glob('*.sh')

        fileList.sort(cmp=None, key=None, reverse=False)
        os.system("mkdir " + logDir)

        i = 0
        while True:
                temp = os.path.join(folder, "temp")
                os.system("qstat > " + temp)
                fileLine = open(os.path.join(folder, temp))
                lineCNT = len(fileLine.readlines())

                if i == len(fileList):
                        print "qsub end..."
                        break
                elif int(lineCNT) > 1000:
                        print "Cluster is full... wait"
                        time.sleep(60)
                else:
                        #print fileList[i]
                        os.system("qsub -V -o " + logDir + " -e " + logDir + ' ' + fileList[i])
                        i = i + 1

#logDir = sys.argv[1]
logDir = "/jdata/scratch/myhwang/Cowork/HHWON/SCRIPTs/Log/"
shRUN(logDir)


