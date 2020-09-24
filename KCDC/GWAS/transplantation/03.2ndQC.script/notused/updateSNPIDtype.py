#!/usr/bin/env python
# coding: utf-8

# In[145]:





wdir = "/ADATA/smkim/JG/03.QC_2nd/CASE/"


# In[151]:


bim = open(wdir+"JG.2nd.QC_snpolisher_rmunrelated_rmaffy_rmdup_convert_indel_flip.bim","r")
print("read bim...")
ref = open(wdir+"../INPUTs/Axiom_KOR.annot.extract.addINDEL.Final.REF.txt","r")
print("read.ref...")
out = open(wdir+"updateSNPID.txt","w")
print("open write...")
trash = open(wdir+"trash.txt","w")


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
    col =  line.split('\t')
    key = col[1]
    A1 = col[4]
    A2 = col[5]
    if key in ref_dic:
        if ref_dic[key] == A1:
            out.write(col[1]+"\t"+col[0]+":"+col[3]+"_"+col[4]+"/"+col[5]+"\n")
        elif ref_dic[key] == A2:
            out.write(col[1]+"\t"+col[0]+":"+col[3]+"_"+col[5]+"/"+col[4]+"\n")
        else:
            trash.write(line + "\n")
bim.close()
out.close()
trash.close()
ref.close()


# In[150]:


print("bye")


# In[ ]:





# In[ ]:




