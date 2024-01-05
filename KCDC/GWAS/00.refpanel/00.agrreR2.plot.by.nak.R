library(ggplot2)
library(dplyr)
library(data.table)
library(iotools)
library(Cairo)
library(cowplot)
library(ggridges)
library(karyoploteR)
library(UpSetR)
library(ggVennDiagram)
library(ggExtra)
library(plotly)
library(ggpubr)
library(ggpp)
#
setwd("/media/epigenome/T7/2023_work_ref/krg2/progress4/13/fig")
#setwd("E:/2023_work_ref/krg2/progress4/13/fig")
#
pop <- "KOR"
#
makeResults <- function(pop){
print(pop)
#
krg <- read.delim(file=paste0("mimicKBA_",pop,"_KRG_GRCh37_AggreR2.RSquare"), head=F, sep="\t")
  colnames(krg) <- c("ID","AF","nSamples","ImputationR2","ValidationAF","ImputationAF")
china <- read.delim(file=paste0("mimicKBA_",pop,"_ChinaMAP_GRCh37_AggreR2.RSquare"), head=F, sep="\t")
  colnames(china) <- c("ID","AF","nSamples","ImputationR2","ValidationAF","ImputationAF")
nard <- read.delim(file=paste0("mimicKBA_",pop,"_NARD_GRCh37_AggreR2.RSquare"), head=F, sep="\t")
  colnames(nard) <- c("ID","AF","nSamples","ImputationR2","ValidationAF","ImputationAF")
gasp <- read.delim(file=paste0("mimicKBA_",pop,"_GAsP_GRCh37_AggreR2.RSquare"), head=F, sep="\t")
  colnames(gasp) <- c("ID","AF","nSamples","ImputationR2","ValidationAF","ImputationAF")
hrc <- read.delim(file=paste0("mimicKBA_",pop,"_HRC_GRCh37_AggreR2.RSquare"), head=F, sep="\t")
  colnames(hrc) <- c("ID","AF","nSamples","ImputationR2","ValidationAF","ImputationAF")
topmed <- read.delim(file=paste0("mimicKBA_",pop,"_TOPMed_GRCh37_AggreR2.RSquare"), head=F, sep="\t")
  colnames(topmed) <- c("ID","AF","nSamples","ImputationR2","ValidationAF","ImputationAF")
#
wellImp <- function(dat){
  print("Identifying well-imputed variatns...")
  dat$wellImputed <- NA
  dat$wellImputed <- ifelse(dat$ImputationR2 >= 0.8, "WellImputed", "Imputed")
  print("Done...")
  return(dat)
}
#
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
  dat2 <- dat %>% group_by(Group,AFBin,wellImputed) %>% summarise(Cnt=length(ID), R2Mean=mean(ImputationR2), R2SD=sd(ImputationR2))
  print("Done")
  return(dat2)
}
#
R2Stat2 <- function(dat){
  print("Calculating Stats for AggreR2...")
  dat2 <- dat %>% group_by(Group,AFBin) %>% summarise(Cnt=length(ID), R2Mean=mean(ImputationR2), R2SD=sd(ImputationR2))
  print("Done")
  return(dat2)
}
#
R2Stat3 <- function(dat){
  print("Calculating Stats for WellImputed...")
  total <- dat %>% group_by(Group,AFBin2) %>% summarise(Cnt=length(ID), R2Mean=mean(ImputationR2), R2SD=sd(ImputationR2)) %>% mutate(wellImputed="Total", Portion=100)
    total <- total[,c("Group","AFBin2","wellImputed","Cnt","R2Mean","R2SD","Portion")]
  imp <- dat %>% filter(wellImputed=="Imputed") %>% group_by(Group,AFBin2,wellImputed) %>% summarise(Cnt=length(ID), R2Mean=mean(ImputationR2), R2SD=sd(ImputationR2))
    imp$Portion <- round((imp$Cnt/total$Cnt)*100, digits=2)
  well <- dat %>% filter(wellImputed=="WellImputed") %>% group_by(Group,AFBin2,wellImputed) %>% summarise(Cnt=length(ID), R2Mean=mean(ImputationR2), R2SD=sd(ImputationR2))
    well$Portion <- round((well$Cnt/total$Cnt)*100, digits=2)
  dat2 <- rbindlist(list(imp,well,total))
  print("Done")
  return(dat2)
}
#
krg.temp <- krg %>% wellImp() %>% AFBinR2(.,"ValidationAF") %>% AFBinWell(.,"ValidationAF") %>% mutate(Group="KRG")
china.temp <- china %>% wellImp() %>% AFBinR2(.,"ValidationAF") %>% AFBinWell(.,"ValidationAF") %>% mutate(Group="ChinaMAP")
nard.temp <- nard %>% wellImp() %>% AFBinR2(.,"ValidationAF") %>% AFBinWell(.,"ValidationAF") %>% mutate(Group="NARD")
gasp.temp <- gasp %>% wellImp() %>% AFBinR2(.,"ValidationAF") %>% AFBinWell(.,"ValidationAF") %>% mutate(Group="GAsP")
#hrc.temp <- hrc %>% wellImp() %>% AFBinR2(.,"ValidationAF") %>% AFBinWell(.,"ValidationAF") %>% mutate(Group="HRC")
topmed.temp <- topmed %>% wellImp() %>% AFBinR2(.,"ValidationAF") %>% AFBinWell(.,"ValidationAF") %>% mutate(Group="TOPMed")
head(krg.temp)

#write.table(krg.temp, file=paste0("imputation_stats_",pop,"_krg.tab"), row.names=F, col.names=T, sep="\t", quote=F)
#write.table(china.temp, file=paste0("imputation_stats_",pop,"_china.tab"), row.names=F, col.names=T, sep="\t", quote=F)
#write.table(nard.temp, file=paste0("imputation_stats_",pop,"_nard.tab"), row.names=F, col.names=T, sep="\t", quote=F)
#write.table(gasp.temp, file=paste0("imputation_stats_",pop,"_gasp.tab"), row.names=F, col.names=T, sep="\t", quote=F)
#write.table(hrc.temp, file=paste0("imputation_stats_",pop,"_hrc.tab"), row.names=F, col.names=T, sep="\t", quote=F)
#write.table(topmed.temp, file=paste0("imputation_stats_",pop,"_topmed.tab"), row.names=F, col.names=T, sep="\t", quote=F)
#
krg.temp2 <- krg.temp %>% R2Stat1()
krg.temp3 <- krg.temp %>% R2Stat2()
krg.temp4 <- krg.temp %>% R2Stat3()
china.temp2 <- china.temp %>% R2Stat1()
china.temp3 <- china.temp %>% R2Stat2()
china.temp4 <- china.temp %>% R2Stat3()
nard.temp2 <- nard.temp %>% R2Stat1()
nard.temp3 <- nard.temp %>% R2Stat2()
nard.temp4 <- nard.temp %>% R2Stat3()
gasp.temp2 <- gasp.temp %>% R2Stat1()
gasp.temp3 <- gasp.temp %>% R2Stat2()
gasp.temp4 <- gasp.temp %>% R2Stat3()
#hrc.temp2 <- hrc.temp %>% R2Stat1()
#hrc.temp3 <- hrc.temp %>% R2Stat2()
#hrc.temp4 <- hrc.temp %>% R2Stat3()
topmed.temp2 <- topmed.temp %>% R2Stat1()
topmed.temp3 <- topmed.temp %>% R2Stat2()
topmed.temp4 <- topmed.temp %>% R2Stat3()
#
merge1 <- rbindlist(list(krg.temp2,china.temp2,nard.temp2,gasp.temp2,hrc.temp2,topmed.temp2))
  merge1$Group <- factor(merge1$Group, levels=c("KRG","NARD","ChinaMAP","GAsP","HRC","TOPMed"))
merge2 <- rbindlist(list(krg.temp3,china.temp3,nard.temp3,gasp.temp3,hrc.temp3,topmed.temp3))
  merge2$Group <- factor(merge2$Group, levels=c("KRG","NARD","ChinaMAP","GAsP","HRC","TOPMed"))
merge3 <- rbindlist(list(krg.temp4,china.temp4,nard.temp4,gasp.temp4,hrc.temp4,topmed.temp4))
  merge3$Group <- factor(merge3$Group, levels=c("KRG","NARD","ChinaMAP","GAsP","HRC","TOPMed"))
head(merge1)
head(merge2)
head(merge3)
write.table(merge1, file=paste0("01.imputation_stats_",pop,"_merge1.tab"), row.names=F, col.names=T, sep="\t", quote=F)
write.table(merge2, file=paste0("01.imputation_stats_",pop,"_merge2.tab"), row.names=F, col.names=T, sep="\t", quote=F)
write.table(merge3, file=paste0("01.imputation_stats_",pop,"_merge3.tab"), row.names=F, col.names=T, sep="\t", quote=F)
#

pop <- "KOR"
#pop <- "EUR"
#pop <- "EAS"
#pop <- "AMR"
#pop <- "AFR"
merge1 <- read.table(paste0("01.imputation_stats_",pop,"_merge1.tab"), header=T, sep="\t") %>% filter(Group!="HRC")
merge2 <- read.table(paste0("01.imputation_stats_",pop,"_merge2.tab"), header=T, sep="\t") %>% filter(Group!="HRC")
merge3 <- read.table(paste0("01.imputation_stats_",pop,"_merge3.tab"), header=T, sep="\t") %>% filter(Group!="HRC")

#

CairoPNG(filename=paste0("01.mimicKBA_",pop,"_aggreR2_plot1_total.png"), width = 750, height = 600, bg="white", pointsize=18)
gg <- merge2
#  gg <- rbind(gg, data.frame(Group=c("KRG","NARD","ChinaMAP","GAsP","HRC","TOPMed"),AFBin=rep("1",6),Cnt=NA,R2Mean=NA,R2SD=NA))
  gg <- rbind(gg, data.frame(Group=c("KRG","NARD","ChinaMAP","GAsP","TOPMed"),AFBin=rep("1",6),Cnt=NA,R2Mean=NA,R2SD=NA))
p <- ggplot(data=gg, aes(x=AFBin, y=R2Mean, group=Group, color=Group, fill=Group)) + geom_line(position=position_nudge(x=0.5,y=0), linetype="solid", size=2) +
        geom_point(position=position_nudge(x=0.5, y=0), aes(shape=Group), size=4) +
        geom_linerange(aes(ymin=R2Mean-R2SD,ymax=R2Mean+R2SD), linetype=3, size=1.5, position=ggpp::position_dodgenudge(x=0.5,y=0,width=0.7,direction="split.x")) +
            ylab(label=bquote(bold("Aggregated R"^2))) + xlab(label="Non-reference allele frequency (%)") +
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
print(p)
dev.off()
#
CairoPNG(filename=paste0("01.mimicKBA_",pop,"_aggreR2_plot1_well.png"), width = 750, height = 600, bg="white", pointsize=18)
gg <- merge1 %>% filter(wellImputed=="WellImputed")
#  gg <- rbind(gg, data.frame(Group=c("KRG","NARD","ChinaMAP","GAsP","HRC","TOPMed"),AFBin=rep("1",6),wellImputed=NA,Cnt=NA,R2Mean=NA,R2SD=NA))
  gg <- rbind(gg, data.frame(Group=c("KRG","NARD","ChinaMAP","GAsP","TOPMed"),AFBin=rep("1",6),wellImputed=NA,Cnt=NA,R2Mean=NA,R2SD=NA))
p <- ggplot(data=gg, aes(x=AFBin, y=R2Mean, group=Group, color=Group, fill=Group)) + geom_line(position=position_nudge(x=0.5,y=0), linetype="solid", size=2) +
    geom_point(position=position_nudge(x=0.5, y=0), aes(shape=Group), size=4) +
    geom_linerange(aes(ymin=R2Mean-R2SD,ymax=R2Mean+R2SD), linetype=3, size=1.5, position=ggpp::position_dodgenudge(x=0.5,y=0,width=0.7,direction="split.x")) +
    ylab(label=bquote(bold("Aggregated R"^2))) + xlab(label="Non-reference allele frequency (%)") +
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
  print(p)
dev.off()
#
CairoPNG(filename=paste0("01.mimicKBA_",pop,"_aggreR2_plot1_noWell.png"), width = 750, height = 600, bg="white", pointsize=18)
gg <- merge1 %>% filter(wellImputed=="Imputed")
#  gg <- rbind(gg, data.frame(Group=c("KRG","NARD","ChinaMAP","GAsP","HRC","TOPMed"),AFBin=rep("1",6),wellImputed=NA,Cnt=NA,R2Mean=NA,R2SD=NA))
gg <- rbind(gg, data.frame(Group=c("KRG","NARD","ChinaMAP","GAsP","TOPMed"),AFBin=rep("1",6),wellImputed=NA,Cnt=NA,R2Mean=NA,R2SD=NA))
p <- ggplot(data=gg, aes(x=AFBin, y=R2Mean, group=Group, color=Group, fill=Group)) + geom_line(position=position_nudge(x=0.5,y=0), linetype="solid", size=2) +
    geom_point(position=position_nudge(x=0.5, y=0), aes(shape=Group), size=4) +
    geom_linerange(aes(ymin=R2Mean-R2SD,ymax=R2Mean+R2SD), linetype=3, size=1.5, position=ggpp::position_dodgenudge(x=0.5,y=0,width=0.7,direction="split.x")) +
    ylab(label=bquote(bold("Aggregated R"^2))) + xlab(label="Non-reference allele frequency (%)") +
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
  print(p)
dev.off()
#
CairoPNG(filename=paste0("01.mimicKBA_",pop,"_wellImputed_plot1_.png"), width = 750, height = 600, bg="white", pointsize=18)
gg <- merge3 %>% filter(wellImputed=="WellImputed")
p <- ggplot(data=gg, aes(x=AFBin2, y=Portion, group=Group, color=Group, fill=Group)) + geom_line(linetype="solid", size=2) +
  geom_point(aes(shape=Group), size=4) +
  ylab(label="% of  well-imputed variants") + xlab(label="Non-reference allele frequency (%)") +
  scale_x_discrete(expand=c(0.05,0),labels=c("AF<0.2","0.2≤AF<0.5","0.5≤AF<5","AF≥5")) +
  scale_y_continuous(breaks=seq(0,100,by=10), limits=c(0,100)) +
  theme(axis.ticks.x.bottom=element_line(size=0),
        panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(colour="white"),
        axis.line=element_line(size=1, colour="#7C7C7C", linetype="solid"),
        axis.title=element_text(colour="Black", face="bold", size=24),
        axis.text=element_text(colour="Black", face="bold", size=24),
        legend.title=element_blank(),
        legend.text=element_text(colour="Black", face="bold", size=20),
        legend.key=element_rect(fill="white"),
        legend.position=c(0.8,0.15))
print(p)
# ggsave(temp_plot, file=paste0("plot_", i,".png"), width = 14, height = 10, units = "cm")
dev.off()


}

#
works <- list("EUR","EAS","AMR","AFR")
  lapply(works, makeResults)