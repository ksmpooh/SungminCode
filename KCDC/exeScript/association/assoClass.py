'''
Created on 2019. 5. 31.

@author: myhwang


Code study..smkim

'''

import os, glob, gzip

from fileInOut import fileInOut  #from fileInOut, load fileInOut class 
fileIO = fileInOut()

class assoClass(object):
	print "CLASS: assoClass"
		
	
	def qtRUN(self, shDir, filDir, traitDir, pedIn, covType, phenoType, epactsTool, chipType):
		print "CLASS: Running qtRUN()..."
		
		for i in range(1, 23):
            
            #get vcfData file in filDir,
            vcfData = glob.glob(filDir + "chr" + str(i) + "_*" + chipType + "*.vcf.gz.tbi")
			#ex)chr1_***_v1****.vcf.gz.tbi
            #chr1_123123123_v1231231231231.vcf.gz.tbi
			
            
            vcfList = [r.replace('\r', '').replace('\n', '') for r in vcfData]
			
			for j in vcfList:
				vcfIn = j.replace(".tbi", '')
                
                #필요한 file name만 남긴다.. _ 로 split하고 chr3,시작번호,마지막번호로 나눈다
				regionArray = vcfIn.replace(filDir, '').replace("_V1_annoINFO_filINFO0.8.vcf.gz", '').replace("_V2_annoINFO_filINFO0.8.vcf.gz", '').split('_')
                
                #chr number와 시작,마지막번호만 남는다.
                #ex) 3:1000-3000
				region = regionArray[0].replace("chr", '') + ':' + regionArray[1] + '-' + regionArray[2]
                
				runType = "_q.linear_" + phenoType
						
				assoOut = vcfIn.replace(filDir, traitDir).replace("_annoINFO_filINFO0.8.vcf.gz",  runType)
				shOut = vcfIn.replace(filDir, shDir).replace("_annoINFO_filINFO0.8.vcf.gz",  runType + "_assoEPACTs.sh")
                
                #make shell script to use epacts
				with open(shOut, 'w') as shWrite:
					shWrite.write(epactsTool + " single --vcf " + vcfIn + 
								" --ped " + pedIn + " --pheno " + phenoType + 
								" --test q.linear --run 8 --field DS --min-mac 5 -min-callrate 0.95 -no-plot" + 
								" --missing NA --out " + assoOut + " --region " + region + covType + '\n')
	
	
	def assoMERGE(self, traitDir, mergeDir, chunkIn, phenoType):
		print "CLASS: Running assoMERGE()..."
		
		chunkData = fileIO.fileRead(chunkIn)
		for i in range(1, 3):
			chipType = "V" + str(i)
			mergeOut = mergeDir + phenoType + '_' + chipType + ".txt"
			with open(mergeOut, 'w') as mergeWrite:
				for j in chunkData[1:]:
					chunkArray = j.split()
					epactsIn = (traitDir + chunkArray[1] + '_' + chipType + "_q.linear_" + phenoType + ".epacts.gz")
					
					if os.path.isfile(epactsIn):
						mergeData = fileIO.gzRead(epactsIn)
						if j == chunkData[1]:
							mergeWrite.write('\n'.join(mergeData[0:]) + '\n')
						else:
							mergeWrite.write('\n'.join(mergeData[0:]) + '\n')
					else:
						os.system("rm -rf " + mergeOut)
						break
						
		