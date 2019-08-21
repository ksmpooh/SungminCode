#!/usr/bin/env python
# coding: utf-8

# In[116]:


'''
BMI association 20190820

KCHIP 136K

DM noDM for V1, V2
'''


# In[145]:


import glob, os


# In[143]:


#function 
def fileRead(fileIn):
    fileIO = open(fileIn,'r')
    InData = [f.replace('/r','').replace('\n',"").replace('\t',' ') for f in fileIO]
    
    fileIO.close()
    return InData

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


# In[138]:



            


# In[4]:


def make_assoSh(version,pedIn,chunksSplit,traitDir,shDir,phenotype,concept,vcfDir):
    for region_list in chunksSplit:
        region = region_list[0].replace("chr","")+":"+region_list[1]+"-"+region_list[2]
        vcfData = vcfDir + region_list[0]+ "_"+region_list[1] +"_"+region_list[2]+"_"+version+"_annoINFO_fileINFO0.8.vcf.gz"
        runType = "_q.linear_" + phenotype+"_"+concept
        
        assoOut = vcfData.replace(vcfDir,traitDir).replace("_annoINFO_fileINFO0.8.vcf.gz",runType)
        shOut = vcfData.replace(vcfDir,shDir).replace("_annoINFO_fileINFO0.8.vcf.gz",runType + "_assoEPACTs.sh")

        
        with open(shOut, 'w') as shWrite:
            shWrite.write("epacts single --vcf "+ vcfData + " --ped " + pedIn +
                          " --pheno "+ phenotype + " --test q.linear --run 8 --field DS --min-mac 5 -min-callrate 0.95 -no-plot"+
                          " --missing NA --out "+ assoOut+" --region " + region + " \n")


# In[2]:


def main():
    wdir = "~/DATA/smkim/KCHIP_130K/BMI_asso/"
    vcfDir = "/LaCie/ghyoon/OAS/"
    outDir = wdir+"RESULTs/"
    inDir = wdir+"INPUTs/"
    scriptDir = wdir+ "SCRIPTs/"
    #vcf file = chr1_999999_999999_V1_annoINFO_fileINFO0.8.vcf.gz
    DM = "KCHIP130K_BMI_adj_DM_20190820"
    noDM = "KCHIP130K_BMI_adj_noDM_20190820"
    Chunk = inDir+"imputation.IMPUTE4.POS.50K_20181114_Final.txt"
    
    
    
    versions = ["V1","V2"]
    concept = "DM"
    
    phenotype = "bmi_inv"
    
    
    Ori_chunksfile = [f.split(" ") for f in fileRead(wdir + Chunk)]
    chunks = [r[1] for r in Ori_chunksfile[1:]]
    chunksSplit = [r.split("_") for r in chunks]
    
    
    traitDir = outDir+"assoRESULTs/"+concept+phenotype+"/"
    make_dir(traitDir)
    shDir = scriptDir+"assoQT_"+phenotype+"/"
    make_dir(shDir)

    
    
    pedfile = "KCHIP130K_BMI_adj_"+concept+"_20190820"
    pedIn = indir + pedfile
    make_assoSh(version="V1",pedIn,chunkSplit,traitDir,shDir,phenotype,concept,vcfDir)


# In[ ]:





# In[ ]:




