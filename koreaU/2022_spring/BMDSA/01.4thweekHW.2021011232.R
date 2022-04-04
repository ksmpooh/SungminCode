#install.packages('rio')

library(tidyverse)
d = rio::import('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6015533/bin/NIHMS957592-supplement-1.xlsx')
class(d)
head(d)
head(d$PatientID)
class(d$TrueRecurrence)
summary(d)
m = as.matrix(d)
head(m[,'TrueRecurrence'])
class(d$SeizureOnsetDays)
class(d$SeizureOnsetDays)
head(d$Seizures)
d1 = d[d$Seizures=='Y',]

head(d1$SeizureOnsetDays)
d1$SeizureOnsetDays2 <- as.numeric(d1$SeizureOnsetDays)
head(is.na(d1$SeizureOnsetDays2))
d1$SeizureOnsetDays[13]
d1[is.na(d1$SeizureOnsetDays2),]$SeizureOnsetDays
# gsub('pattern in your character', 'new character you want to replace', vectors for your character)
d1$SeizureOnsetDays3 <- gsub('<', '', d1$SeizureOnsetDays)
head(d1$SeizureOnsetDays3)
min(d1$SeizureOnsetDays3)

head(d$p.Protein)
table(d$p.Protein)

gsub('[^0-9]', '', 'p.R102X')
for (i in 1:293) {
  a <- gsub('[^0-9]', '', d[i, 'p.Protein'] )
}
a
head(d1)
gsub('[^0-9]', '', d$p.Protein)


#3-1
d$p.Protein_pos <- gsub('[^0-9]', '', d$p.Protein)
#3-2
c <- as.data.frame(table(d$p.Protein))
c
## we can check the number of unique variants in the dataset by:
c %>% filter(Freq == 1) %>% count() # 163 unique variants
c %>% filter(Freq != 1) %>% count() # 40 multiple variants
##How many unique variants you can find? and which variants are occurred in multiple times?
#163 unique variants, 40 variants in multiple times

c1 <- c %>% filter(Freq != 1)
## please check the recurrent mutations for each group.
d %>% filter(p.Protein %in% c1$Var1) %>% 
  group_by(p.Protein,Classification) %>%
  summarise(n=n()) %>% filter(n == 1)


table(d$TrueRecurrence)
#3-3
## 1) calculate the proportaion of the phenotypes
table(d$Classification)/nrow(d)

##  2) to whether females and males have different occurrence in each disorder
### other value change to "unkwnown" except F/M in PatientSex Column
d %>% mutate(New_PatientSex = ifelse(PatientSex == "F" | PatientSex == "M",PatientSex,"unknown")) %>% 
  group_by(Classification,New_PatientSex) %>%
  summarise(n = n()) %>%
  mutate(freq=n/sum(n))

## 3) Let's find out how many mutation consequences are observed in each phenotype. 
d %>% group_by(Classification,Effect) %>%
  summarise(n=n())

#3-4
table(d$Effect)
#If you checked the recurrent mutations, you might want to find a locus where two or more variants occur. 
#Such loci might indicate functionally important position of the gene and you might find some insight as to a cause of disease.
d1 <- d %>% group_by(c.DNA,Effect) %>%
  summarise(n=n()) %>% filter(n != 1) %>%
  mutate(DNA_pos = as.integer(gsub('[^0-9]', '', c.DNA))) %>%
  arrange(DNA_pos) #%>%

table(d1$Effect)

d1 %>%  ggplot(aes(Effect)) + 
  geom_bar(position = fill)

d %>% group_by(p.Protein,Effect) %>%
  summarise(n=n()) %>% filter(n != 1) %>%
  mutate(Protein_pos = as.integer(gsub('[^0-9]', '', p.Protein))) %>%
  arrange(Protein_pos)

d %>% select(c.DNA,p.Protein,p.Protein_pos,Effect) %>%
  filter(p.Protein %in% c1$Var1) %>% mutate(p.Protein_pos = as.integer(p.Protein_pos)) %>%
  arrange(p.Protein_pos)

d %>% select(c.DNA,p.Protein,p.Protein_pos,Effect) %>%
  filter(p.Protein %in% c1$Var1) %>% 
  group_by(c.DNA,p.Protein) %>%
  summarise(n=n()) %>% 
  mutate(Protein_pos = as.integer(gsub('[^0-9]', '', p.Protein))) %>%
  arrange(Protein_pos)

# 3-5
d %>% select(Effect,Classification) %>%
  ggplot(aes(x=Classification,fill=Effect)) +
  geom_bar(position='fill')
