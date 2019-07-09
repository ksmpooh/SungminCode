'''
Created on 2017. 6. 1.
@author: HMY

Code study..smkim on 2019.07.09
make shell script using python
'''

#basic library in python
import os, sys, glob, time



def shRUN(logDir):
        folder = os.getcwd() # current working directory
        fileList = glob.glob('*.sh') # all .sh file put in filelist 

        fileList.sort(cmp=None, key=None, reverse=False) #sort
        os.system("mkdir " + logDir) #mkdir log directory

        i = 0
        while True:
                temp = os.path.join(folder, "temp") #temp file을 만들고
                os.system("qstat > " + temp) #qstat(working queue 확인) -> temp
                fileLine = open(os.path.join(folder, temp))
                lineCNT = len(fileLine.readlines()) #line count on fileline

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


