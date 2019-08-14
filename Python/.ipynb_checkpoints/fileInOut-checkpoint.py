#file IO class
# in linux command
import os,gzip

class fileInOut():
    def makeDir(self, outDir):
        print "makeDir()"
        
        os.system("mkdir " + outDir)
    
    
    def fileRead(self,fileIn):
        fileList = open(fileIn,"r")
        inData = [f.replace('')]