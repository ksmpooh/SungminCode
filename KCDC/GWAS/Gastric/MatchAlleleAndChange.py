#!/usr/bin/env python
# coding: utf-8

# In[145]:


wdir = "c:/Users/user/Desktop/KCDC/Gastric/Ref/"


# In[151]:


bim = open(wdir+"CASE_CONTROL_MERGE_rmfreq.bim","r")
ref = open(wdir+"Axiom_KOR.annot.extract.addINDEL.Final.REF.txt","r")
out = open(wdir+"new_Case_Control.bim","w")
trash = open(wdir+"trash.txt","w")


# In[152]:


ref_dic = {}


# In[153]:


while True:
    line = ref.readline().replace("\n","")
    if not line : break
    col  = line.split("\t")
    ref_dic[col[0]] = col[1]


# In[154]:


while True:
    line = bim.readline().replace("\n","")
    if not line : break
    col =  line.split('\t')
    key = col[1]
    A1 = col[4]
    A2 = col[5]
    if key in ref_dic:
        if ref_dic[key] == A1:
            out.write(col[0]+"\t"+col[0]+":"+col[3]+"_"+col[4]+"/"+col[5]+"\t"+col[2]+"\t"+col[3]+"\t"+col[4]+"\t"+col[5]+"\n")
        elif ref_dic[key] == A2:
            out.write(col[0]+"\t"+col[0]+":"+col[3]+"_"+col[5]+"/"+col[4]+"\t"+col[2]+"\t"+col[3]+"\t"+col[5]+"\t"+col[4]+"\n")
        else:
            trash.write(line + "\n")
bim.close()
out.close()
trash.close()
ref.close()


# In[150]:


print("hi")


# In[ ]:





# In[ ]:




