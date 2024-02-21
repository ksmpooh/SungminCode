## 8K imputation
library(tidyverse)
library(stringr)
library(ggplot2)

library(cowplot)
library(UpSetR)

library(ggVennDiagram)
library(ggExtra)
library(plotly)
library(ggpubr)
library(ggpp)

#library(karyoploteR)
#library(ggridges)
library(Cairo)


setwd("~/Desktop/KCDC/imputation/8K/INFO/")

flist = grep(list.files("./"),pattern = ".info", value=TRUE)
flist

out <- NULL

for (i in 1:length(flist)) {
  df <- read_table(flist[i])
  df <- df %>% select(SNP,MAF,Rsq,Genotyped) %>% mutate(chrom=str_split_fixed(SNP,":",3)[,1])
  out <- rbind(out,df)
}

#><
#ori <- out
out <- ori
head(out)
out %>% filter(MAF >= 0.05) %>% head()
out %>% filter(Genotyped == "Imputed") %>% select(-Genotyped) %>%
  mutate(MAF = ifelse(MAF >= 0.05,"MAF>5",ifelse(MAF>=0.01,"1<MAF<5",ifelse(MAF >=0.005,"0.5<MAF<1"))))

#out %>% filter(Genotyped == "Imputed") %>% select(-Genotyped) -> out
head(out)

#out %>% filter(Genotyped == "Imputed") %>% select(-Genotyped) %>% write.table("KCHIP130K_IMP8K_allCHR_variantINFO.txt",col.names = T,row.names = F,quote = F,sep = "\t")


wellImp <- function(dat){
  print("Identifying well-imputed variatns...")
  dat$wellImputed <- NA
  dat$wellImputed <- ifelse(dat$Rsq >= 0.8, "WellImputed", "Imputed")
  print("Done...")
  return(dat)
}

AFBinR2 <- function(dat, col){
  print("Identifying AF bin for R2...")
  print(col)
  idx <- which(colnames(dat) == col)
  af <- dat[ ,idx]
  dat$AFBin <- NA
  dat$AFBin <- ifelse( af < 0.002, "<0.002", dat$AFBin)
  dat$AFBin <- ifelse( 0.002 <= af & af < 0.005, "0.002~0.005", dat$AFBin)
  dat$AFBin <- ifelse( 0.005 <= af & af < 0.01, "0.005~0.01", dat$AFBin)
  dat$AFBin <- ifelse( 0.01 <= af & af < 0.02, "0.01~0.02", dat$AFBin)
  dat$AFBin <- ifelse( 0.02 <= af & af < 0.05, "0.02~0.05", dat$AFBin)
  dat$AFBin <- ifelse( 0.05 <= af & af < 0.1, "0.05~0.1", dat$AFBin)
  dat$AFBin <- ifelse( 0.1 <= af & af < 0.2, "0.1~0.2", dat$AFBin)
  dat$AFBin <- ifelse( 0.2 <= af & af < 0.3, "0.2~0.3", dat$AFBin)
  dat$AFBin <- ifelse( 0.3 <= af & af < 0.4, "0.3~0.4", dat$AFBin)
  dat$AFBin <- ifelse( 0.4 <= af & af < 0.5, "0.4~0.5", dat$AFBin)
  dat$AFBin <- ifelse( 0.5 <= af, "0.5<=", dat$AFBin)
  temp <- c("<0.002","0.002~0.005","0.005~0.01","0.01~0.02","0.02~0.05","0.05~0.1","0.1~0.2","0.2~0.3","0.3~0.4","0.4~0.5","0.5<=")
  dat$AFBin <- factor(dat$AFBin, levels=temp)
  print("Done...")
  return(dat)
}

AFBinR2_v2 <- function(dat, col){
  print("Identifying AF bin for R2...")
  print(col)
  idx <- which(colnames(dat) == col)
  af <- dat[ ,idx]
  dat$AFBin <- NA
  dat$AFBin <- ifelse( af < 0.002, "<0.002", dat$AFBin)
  dat$AFBin <- ifelse( 0.002 <= af & af < 0.005, "0.002~0.005", dat$AFBin)
  dat$AFBin <- ifelse( 0.005 <= af & af < 0.01, "0.005~0.01", dat$AFBin)
  dat$AFBin <- ifelse( 0.01 <= af & af < 0.02, "0.01~0.02", dat$AFBin)
  dat$AFBin <- ifelse( 0.02 <= af & af < 0.05, "0.02~0.05", dat$AFBin)
  dat$AFBin <- ifelse( 0.05 <= af & af < 0.1, "0.05~0.1", dat$AFBin)
  #dat$AFBin <- ifelse( 0.1 <= af & af < 0.2, "0.1~0.2", dat$AFBin)
  dat$AFBin <- ifelse( 0.1 <= af, "0.1<=", dat$AFBin)
  temp <- c("<0.002","0.002~0.005","0.005~0.01","0.01~0.02","0.02~0.05","0.05~0.1","0.1<=")
  dat$AFBin <- factor(dat$AFBin, levels=temp)
  print("Done...")
  return(dat)
}

#
AFBinWell <- function(dat, col){
  print("Identifying AF bin for Well...")
  print(col)
  idx <- which(colnames(dat) == col)
  af <- dat[ ,idx]
  dat$AFBin2 <- NA
  dat$AFBin2 <- ifelse( af < 0.002, "<0.002", dat$AFBin2)
  dat$AFBin2 <- ifelse( 0.002 <= af & af < 0.005, "0.002~0.005", dat$AFBin2)
  dat$AFBin2 <- ifelse( 0.005 <= af & af < 0.5, "0.005~0.5", dat$AFBin2)
  dat$AFBin2 <- ifelse( 0.5 <= af, "0.5<=", dat$AFBin2)
  temp <- c("<0.002","0.002~0.005","0.005~0.5","0.5<=")
  dat$AFBin2 <- factor(dat$AFBin2, levels=temp)
  print("Done...")
  return(dat)
}
#
R2Stat1 <- function(dat){
  print("Calculating Stats for AF bin...")
  #dat2 <- dat %>% group_by(Group,AFBin,wellImputed) %>% summarise(Cnt=length(ID), R2Mean=mean(ImputationR2), R2SD=sd(ImputationR2))
  #dat2 <- dat %>% group_by(chrom,AFBin,wellImputed) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  dat2 <- dat %>% group_by(AFBin,wellImputed) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  print("Done")
  return(dat2)
}
#
R2Stat2 <- function(dat){
  print("Calculating Stats for AggreR2...")
  #dat2 <- dat %>% group_by(Group,AFBin) %>% summarise(Cnt=length(ID), R2Mean=mean(ImputationR2), R2SD=sd(ImputationR2))
  #dat2 <- dat %>% group_by(chrom,AFBin) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  dat2 <- dat %>% group_by(AFBin) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  print("Done")
  return(dat2)
}
R2Stat2_chr <- function(dat){
  print("Calculating Stats for AggreR2...chr")
  #dat2 <- dat %>% group_by(Group,AFBin) %>% summarise(Cnt=length(ID), R2Mean=mean(ImputationR2), R2SD=sd(ImputationR2))
  dat2 <- dat %>% group_by(chrom,AFBin) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  #dat2 <- dat %>% group_by(AFBin) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  print("Done")
  return(dat2)
}

#
head(toy)
R2Stat3 <- function(dat){
  print("Calculating Stats for WellImputed...")
  total <- dat %>% group_by(chrom,AFBin2) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq)) %>% mutate(wellImputed="Total", Portion=100)
  total <- total[,c("chrom","AFBin2","wellImputed","Cnt","R2Mean","R2SD","Portion")]
  head(total)
  imp <- dat %>% filter(wellImputed=="Imputed") %>% group_by(chrom,AFBin2,wellImputed) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  imp$Portion <- round((imp$Cnt/total$Cnt)*100, digits=2)
  well <- dat %>% filter(wellImputed=="WellImputed") %>% group_by(chrom,AFBin2,wellImputed) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  well$Portion <- round((well$Cnt/total$Cnt)*100, digits=2)
  dat2 <- rbindlist(list(imp,well,total))
  print("Done")
  return(dat2)
}

R2Stat3 <- function(dat){
  print("Calculating Stats for WellImputed...")
  total <- dat %>% group_by(chrom,AFBin2) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq)) %>% mutate(wellImputed="Total", Portion=100)
  head(total)
  total <- total[,c("chrom","AFBin2","wellImputed","Cnt","R2Mean","R2SD","Portion")]
  head(total)
  imp <- dat %>% filter(wellImputed=="Imputed") %>% group_by(chrom,AFBin2,wellImputed) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  imp$Portion <- round((imp$Cnt/total$Cnt)*100, digits=2)
  well <- dat %>% filter(wellImputed=="WellImputed") %>% group_by(chrom,AFBin2,wellImputed) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq))
  well$Portion <- round((well$Cnt/total$Cnt)*100, digits=2)
  dat2 <- rbindlist(list(imp,well,total))
  print("Done")
  return(dat2)
}
head(df.tmp)
table(df.tmp$AFBin2)
df.tmp %>% group_by(wellImputed,AFBin2) %>% summarise(Cnt=length(SNP), R2Mean=mean(Rsq), R2SD=sd(Rsq)) %>% #head()
  group_by(AFBin2) %>% mutate(Portion = 100*Cnt/sum(Cnt)) -> df.tmp.t3

head(df.tmp.t3)
head(toy)
#out.temp <- out %>% wellImp() %>% AFBinR2(.,"MAF") %>% AFBinWell(.,"MAF")

df <- read_table("KCHIP130K_IMP8K_allCHR_variantINFO.txt")
head(df)
toy <- df %>% filter(chrom==22)
head(toy)

toy %>% wellImp() %>% AFBinR2(.,"MAF") %>% AFBinWell(.,"MAF") -> toy
toy.t1 <- toy %>% R2Stat1()
toy.t2 <- toy %>% R2Stat2()
toy.t3 <- toy %>% R2Stat3()


df %>% wellImp() %>% AFBinR2(.,"MAF") %>% AFBinWell(.,"MAF") -> df.tmp
df.tmp.t1 <- df.tmp %>% R2Stat1()
df.tmp.t2 <- df.tmp %>% R2Stat2()
df.tmp.t2_chr <- df.tmp %>% R2Stat2_chr()
df.tmp.t3 <- df.tmp %>% mutate(chrom = as.factor(chrom)) %>% R2Stat3()
head(df.tmp)
head(df.tmp.t1)
head(df.tmp.t2)
head(df.tmp.t3)

head(df.tmp.t1)
head(df.tmp.t2)
head(df.tmp)

df.tmp.t1 %>% ggplot(aes(x=AFBin,y=R2Mean,group=wellImputed,color=wellImputed,fill=wellImputed)) +
#df.tmp.t1 %>% ggplot(aes(x=AFBin,y=R2Mean)) + 
  geom_line(position=position_nudge(x=0.5,y=0), linetype="solid", size=2) +
  geom_point(position=position_nudge(x=0.5, y=0), size=4) +
  geom_linerange(aes(ymin=R2Mean-R2SD,ymax=R2Mean+R2SD), linetype=3, linewidth=1.5, position=ggpp::position_dodgenudge(x=0.5,y=0,width=0.7,direction="split.x")) +
  #ylab(label=bquote(bold("Aggregated R"^2))) + xlab(label="Non-reference allele frequency (%)") +
  ylab(label=bquote(bold("R"^2))) + xlab(label="Non-reference allele frequency (%)") +
  scale_x_discrete(labels=c("0","0.2","0.5","1","2","5","10","20","30","40","50","100")) +
  scale_y_continuous(expand=c(0, 0), breaks=seq(0,1,by=0.1)) +
  theme(axis.ticks.x.bottom=element_line(size=0),
        panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(colour="white"),
        axis.line=element_line(size=1, colour="#7C7C7C", linetype="solid"),
        axis.title=element_text(colour="Black", face="bold", size=24),
        axis.text=element_text(colour="Black", face="bold", size=24),
        legend.title=element_blank(),
        legend.text=element_text(colour="Black", face="bold", size=20),
        legend.key=element_rect(fill="white"),
        legend.position=c(0.9,0.15))

head(df.tmp.t2_chr)

head(df.tmp)
#df.tmp.t2_chr %>% ggplot(aes(x=AFBin,y=R2Mean,fill=AFBin)) +
df.tmp %>% ggplot(aes(x=AFBin,y=Rsq,fill=AFBin)) + 
  geom_boxplot() + 
  #geom_linerange(aes(ymin=R2Mean-R2SD,ymax=R2Mean+R2SD), linetype=3, size=1.5, position=ggpp::position_dodgenudge(x=0.5,y=0,width=0.7,direction="split.x")) +
  #geom_boxplot(aes(ymin=R2Mean-R2SD,ymax=R2Mean+R2SD), linetype=3, size=1.5, position=ggpp::position_dodgenudge(x=0.5,y=0,width=0.7,direction="split.x")) +
  ylab(label=bquote(bold("R"^2))) + xlab(label="Non-reference allele frequency (%)") +
  scale_x_discrete(labels=c("0","0.2","0.5","1","2","5","10","20","30","40","50","100")) +
  scale_y_continuous(expand=c(0, 0), breaks=seq(0,1,by=0.1)) +
  ylim(c(0.7, 1)) +
  theme(axis.ticks.x.bottom=element_line(size=0),
        #panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(colour="white"),
        axis.line=element_line(size=1, colour="#7C7C7C", linetype="solid"),
        axis.title=element_text(colour="Black", face="bold", size=24),
        axis.text=element_text(colour="Black", face="bold", size=24),
        legend.title=element_blank(),
        legend.text=element_text(colour="Black", face="bold", size=20),
        legend.key=element_rect(fill="white"),
        #legend.position=c(0.9,0.15)
        legend.position="none")







head(df.tmp.t3)
df.tmp.t3 %>% #head()#filter(wellImputed == "WellImputed") %>%
  mutate(wellImputed = ifelse(wellImputed == "WellImputed","R2>0.8",'R2<0.8')) %>%
  ggplot(aes(x=AFBin2, y=Portion, group=wellImputed, color=wellImputed, fill=wellImputed)) +
  #geom_line(linetype="solid", size=2) +
#ggplot(data=gg, aes(x=AFBin2, y=Portion, group=Group, color=Group, fill=Group)) + geom_line(linetype="solid", size=2) +
  geom_bar(position="fill",stat="identity") +
  #geom_text(aes(label = round(Portion/100,2)), position = position_stack(vjust = 0.5))
  geom_text(aes(y = round(Portion/100,2),label=round(Portion/100,2)),size=10,fontface = "bold",color = "white",position = position_stack(vjust = 0.5)) +
#  geom_point( size=4) +
  ylab(label="% of  well-imputed variants") + xlab(label="Non-reference allele frequency (%)") +
  scale_x_discrete(expand=c(0.05,0),labels=c("AF<0.2","0.2<AF<0.5","0.5<AF<5","AF>5")) +
  #scale_y_continuous(breaks=seq(0,100,by=10), limits=c(0,100)) #+
  #scale_y_continuous(breaks=seq(0,1,by=10), limits=c(0,100)) #+
  theme(axis.ticks.x.bottom=element_line(size=0),
        panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(colour="white"),
        axis.line=element_line(size=1, colour="#7C7C7C", linetype="solid"),
        axis.title=element_text(colour="Black", face="bold", size=24),
        axis.text=element_text(colour="Black", face="bold", size=24),
        legend.title=element_blank(),
        legend.text=element_text(colour="Black", face="bold", size=20),
        legend.key=element_rect(fill="white"),
        legend.position=c(0.9,0.15))

df.tmp %>% group_by(chrom,wellImputed) %>% count(chrom) %>% #head()
  pivot_wider(names_from = wellImputed,values_from = n) %>% writexl::write_xlsx("KCHIP130K_IMP8K_wellINFOcount.byCHR.xlsx")


df.tmp %>% group_by(chrom,wellImputed) %>% count(AFBin) %>% #head()
  pivot_wider(names_from = wellImputed,values_from = n) %>% #head()
  mutate(chrom = as.factor(chrom),coverage = WellImputed/(WellImputed+Imputed)) %>% select(-Imputed,-WellImputed) %>% 
  pivot_wider(names_from = AFBin,values_from = coverage) %>% writexl::write_xlsx("KCHIP130K_IMP8K_genomic_overage.byCHR.xlsx")


df %>% wellImp() %>% AFBinR2_v2(.,"MAF") %>% AFBinWell(.,"MAF") %>% group_by(chrom,wellImputed) %>% count(AFBin) %>% #head()
  pivot_wider(names_from = wellImputed,values_from = n) %>% #head()
  mutate(chrom = as.factor(chrom),coverage = WellImputed/(WellImputed+Imputed)) %>% select(-Imputed,-WellImputed) -> df.tmp.t3

head(df.tmp.t3)

df.tmp %>% group_by(wellImputed) %>% count(AFBin) %>% #head()
  pivot_wider(names_from = wellImputed,values_from = n) %>% #head()
  mutate(coverage = WellImputed/(WellImputed+Imputed)) %>% select(-Imputed,-WellImputed) 

df.tmp.t3 %>% filter(AFBin != '0.5<=') %>% head()
df.tmp.t3 %>% #filter(AFBin != '0.5<=') %>%
  ggplot(aes(x=AFBin,y=chrom,fill=coverage)) + 
  geom_tile(width = 0.8, height = 0.8) +
  scale_fill_gradient(low = "red", high = "blue",name = "Genomic\nCoverage") + 
  ylab(label="Chromosome") + xlab(label="Non-reference allele frequency (%)") +
  scale_x_discrete(labels=c("0~0.2","0.2~0.5","0.5~1","1~2","2~5","5~10","10~")) +
  theme(axis.ticks.x.bottom=element_line(size=0),
        panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(colour="white"),
        axis.line=element_line(size=1, colour="#7C7C7C", linetype="solid"),
        axis.title=element_text(colour="Black", face="bold", size=24),
        axis.text=element_text(colour="Black", face="bold", size=20),
        #legend.title=element_blank(),
        legend.title=element_text(face='bold',size=20),
        legend.text=element_text(colour="Black", face="bold", size=20),
        legend.key=element_rect(fill="white"))

