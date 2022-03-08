#fastSIMCOAL

setwd("C:/Users/Bashir/Desktop/Meehania/Population results/FastSimCoal Demographic history/PCA_clusters Fastsimcoal/bestLhood")

all_lineage_geneflow<-read.table("all_lineage_geneflow.lhoods")
ancient_CH1_vs_CH2<-read.table("ancient_CH1_vs_CH2.lhoods")
ancient_CHJP_vs_CH1<-read.table("ancient_CHJP_vs_CH1.lhoods")
ancient_geneflow<-read.table("ancient_geneflow.lhoods")
ancient_vs_CHJP_vs_CH2<-read.table("ancient_vs_CHJP_vs_CH2.lhoods")
no_geneflow<-read.table("no_geneflow.lhoods")
only_recent_geneflow<-read.table("only_recent_geneflow.lhoods")
sub_CH1_vs_CH2<-read.table("sub_CH1_vs_CH2.lhoods")
sub_CHJP_vs_CH2<-read.table("sub_CHJP_vs_CH2.lhoods")
subchinaCHJP_vs_CH1<-read.table("subchinaCHJP_vs_CH1.lhoods")

par(mfrow=c(1,1))

boxplot(range = 0,
        all_lineage_geneflow$V1,
        ancient_CH1_vs_CH2$V1,
        ancient_CHJP_vs_CH1$V1,
        ancient_geneflow$V1,
        ancient_vs_CHJP_vs_CH2$V1,
        no_geneflow$V1,
        only_recent_geneflow$V1,
        sub_CH1_vs_CH2$V1,
        sub_CHJP_vs_CH2$V1,
        subchinaCHJP_vs_CH1$V1,
        ylab="Likelihood",xaxt="n")

lablist.x<-as.vector(c("all_lineage_geneflow",
                       "ancient_CH1_vs_CH2",
                       "ancient_CHJP_vs_CH1",
                       "ancient_geneflow",
                       "ancient_vs_CHJP_vs_CH2",
                       "no_geneflow",
                       "only_recent_geneflow",
                       "sub_CH1_vs_CH2",
                       "sub_CHJP_vs_CH2",
                       "subchinaCHJP_vs_CH1"))
axis(1, at=seq(1, 14, by=1), labels = FALSE)
text(x = seq(1, 14, by=1), par("usr")[3] - 0.2, labels = lablist.x, srt = 45, pos = 1, xpd = TRUE)

library(tidyverse)
#install.packages("ggthemes")
library(ggthemes)

all_lineage_geneflow<-read.table("all_lineage_geneflow.lhoods") %>%
  rename(flow=V1)
all_lineage_geneflow$ID<-rep("all_lineage_geneflow.lhoods",nrow(all_lineage_geneflow))

ancient_CH1_vs_CH2<-read.table("ancient_CH1_vs_CH2.lhoods") %>%
  rename(flow=V1)
ancient_CH1_vs_CH2$ID<-rep("ancient_CH1_vs_CH2.lhoods",nrow(ancient_CH1_vs_CH2))

ancient_CHJP_vs_CH1<-read.table("ancient_CHJP_vs_CH1.lhoods") %>%
  rename(flow=V1)
ancient_CHJP_vs_CH1$ID<-rep("ancient_CHJP_vs_CH1.lhoods",nrow(ancient_CHJP_vs_CH1))

ancient_geneflow<-read.table("ancient_geneflow.lhoods") %>%
  rename(flow=V1)
ancient_geneflow$ID<-rep("ancient_geneflow.lhoods",nrow(ancient_geneflow))

ancient_vs_CHJP_vs_CH2<-read.table("ancient_vs_CHJP_vs_CH2.lhoods") %>%
  rename(flow=V1)
ancient_vs_CHJP_vs_CH2$ID<-rep("ancient_vs_CHJP_vs_CH2.lhoods",nrow(ancient_vs_CHJP_vs_CH2))

no_geneflow<-read.table("no_geneflow.lhoods") %>%
  rename(flow=V1)
no_geneflow$ID<-rep("no_geneflow.lhoods",nrow(no_geneflow))

only_recent_geneflow<-read.table("only_recent_geneflow.lhoods") %>%
  rename(flow=V1)
only_recent_geneflow$ID<-rep("only_recent_geneflow.lhoods",nrow(only_recent_geneflow))

sub_CH1_vs_CH2<-read.table("sub_CH1_vs_CH2.lhoods") %>%
  rename(flow=V1)
sub_CH1_vs_CH2$ID<-rep("sub_CH1_vs_CH2.lhoods",nrow(sub_CH1_vs_CH2))

sub_CHJP_vs_CH2<-read.table("sub_CHJP_vs_CH2.lhoods") %>%
  rename(flow=V1)
sub_CHJP_vs_CH2$ID<-rep("sub_CHJP_vs_CH2.lhoods",nrow(sub_CHJP_vs_CH2))

subchinaCHJP_vs_CH1<-read.table("subchinaCHJP_vs_CH1.lhoods") %>%
  rename(flow=V1)
subchinaCHJP_vs_CH1$ID<-rep("subchinaCHJP_vs_CH1.lhoods",nrow(subchinaCHJP_vs_CH1))

datacombined <- rbind(all_lineage_geneflow,
                      ancient_CH1_vs_CH2,
                      ancient_CHJP_vs_CH1,
                      ancient_geneflow,
                      ancient_vs_CHJP_vs_CH2,
                      no_geneflow,
                      only_recent_geneflow,
                      sub_CH1_vs_CH2,
                      sub_CHJP_vs_CH2,
                      subchinaCHJP_vs_CH1)

qplot( x=ID , y=flow , data=datacombined , geom=c("boxplot","jitter") , colour=ID) +
  theme_minimal()

qplot( x=ID , y=flow , data=datacombined , geom=c("boxplot","jitter") , colour=ID) +
  theme_hc()
