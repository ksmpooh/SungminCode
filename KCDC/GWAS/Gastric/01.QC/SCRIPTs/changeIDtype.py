#!/usr/bin/env python
# coding: utf-8

# In[145]:


wdir = "/DATA/smkim/Gastric/QC_2nd/OUTPUTs/"


# In[151]:


bim = open(wdir+"CASE/KNIH.RAW.Gastric.2nd.rmSNP.rmLQSamples.rmaffy.fil.convert.bim","r")
print("read bim...")
ref = open(wdir+"../INPUTs/Axiom_KOR.annot.extract.addINDEL.Final.REF.txt","r")
print("read.ref...")
out = open(wdir+"../INPUTs/change_SNPID_type_Axiom_KOR.annot.extract.addINDEL.FINAL.REF.txt","w")
print("open write...")
trash = open(wdir+"/nonmatchingID.txt","w")


# In[152]:


ref_dic = {}


# In[153]:


while True:
    line = ref.readline().replace("\n","")
    if not line :
	print("MAke dict...Done")
	break
    col  = line.split("\t")
    ref_dic[col[0]] = col[1]


# In[154]:


while True:
    line = bim.readline().replace("\n","")
    if not line :
	print("New BIM..")
	break
    #print(line)
    col =  line.split('\t')
    key = col[1]
    A1 = col[4]
    A2 = col[5]
    if key in ref_dic:
        if ref_dic[key] == A1:
            out.write(col[0]+":"+col[3]+"_"+col[4]+"/"+col[5]+"\t"+ref_dic[key]+"\n")
        elif ref_dic[key] == A2:
            out.write(col[0]+":"+col[3]+"_"+col[5]+"/"+col[4]+"\t"+ref_dic[key]+"\n")
        else:
            trash.write(line + "\n")
bim.close()
out.close()
trash.close()
ref.close()


# In[150]:


print("bye!")


# In[ ]:





# In[ ]:




