# Requires the "here" library
library(here)
# Requires ggplot2
library(ggplot2)
# A nice ggplot2 theme
source(here("scripts","ggplot_theme.R"))

#--------------------------#
# Load and prepare the data
dataset<-"sp1" # One of "sp1", "sp2", ...
sampl<-"sample1" # One of "sample1", "sample1_2"

# Load the busco results table
rarefaction_dat<-read.table(here(paste("data/",dataset,sep=""),
                                 paste(sampl,"_subsampling_rarefaction_dat.csv",sep="")),
                            header=TRUE,sep=",")

sp_busco_table<-read.table(here("data","sp_busco_table.tab"),header=TRUE,sep="\t")
sp_g_busco<-sp_busco_table$g_busco[which(sp_busco_table$spp==dataset)]

busco_table<-read.table(here(paste("data/",dataset,sep=""),
                             paste(sampl,".polished_busco_full_table.tsv",sep="")),
                        header=FALSE,sep="\t",fill = TRUE)
colnames(busco_table)<-c("busco_id","status","transcript","length","length2")
# get the total nr buscos from the busco_table
busco_total<-length(unique(busco_table$busco_id))
# Subset results to only complete (Single copy) or complete (duplicated) buscos
busco_table<-busco_table[which(busco_table$status=="Complete" | busco_table$status=="Duplicated"),]
busco_score<-(length(unique(busco_table$busco_id))/busco_total)*100

#--------------#
# Plot the data
plot1<-ggplot(data=rarefaction_dat)+
  geom_line(aes(x=reads/100000,y=perc_compl))+
  geom_line(aes(x=reads/100000,y=sd_upr_perc_compl),linetype="dashed")+
  geom_line(aes(x=reads/100000,y=sd_lwr_perc_compl),linetype="dashed")+
  geom_line(aes(x=reads/100000,y=max_perc_compl),linetype="dashed",colour="red")+
  geom_hline(yintercept = busco_score,linetype="dotdash")+ # This line is only relevant for my data, I'm
  geom_label(x=-Inf,y=busco_score,hjust=0,label="IsoSeq3 Clusters")+ # This line is only relevant for my data, I'm
  geom_hline(yintercept = sp_g_busco,linetype="dotdash")+ # This line is only relevant for my data, I'm
  geom_label(x=-Inf,y=sp_g_busco,hjust=0,label="Genome Assembly")+ # This line is only relevant for my data, I'm
  ylim(0,100)+
  #  xlim(0,40)+
  #  geom_point(aes(x=size,y=mean,colour=data))+
  xlab(expression(paste("# Subsampled FL CCS Reads (x 10"^5,")")))+
  ylab("% Complete BUSCO Genes")+
  ggtitle(label=paste(dataset,"-",sampl))+
  my.theme+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))

ggsave(filename=here("figures",paste(dataset,"_plot1.png",sep="")),device = "png",dpi = 300,width = 20,height=15,units="cm")
plot1
dev.off()
ggsave(filename=here("figures",paste(dataset,"_plot1.pdf",sep="")),device = "pdf",dpi = 300,width = 20,height=15,units="cm")
plot1
dev.off()

thresh<-min(tail(rarefaction_dat$reads,n=15))

plot2<-ggplot()+
  geom_point(data=rarefaction_dat[which(rarefaction_dat$reads>thresh),],
                                      aes(x=reads/100000, y=perc_compl))+
  geom_smooth(data=rarefaction_dat[which(rarefaction_dat$reads>thresh),],
              aes(x=reads/100000, y=perc_compl),method = "lm")+
  xlab(expression(paste("# Subsampled FL CCS Reads (x 10"^5,")")))+
  ylab("% Complete BUSCO Genes")+
  ggtitle(label=paste(dataset,"-",sampl))+
  my.theme

ggsave(filename=here("figures",paste(dataset,"_plot2.png",sep="")),device = "png",dpi = 300,width = 10,height=7,units="cm")
plot2
dev.off()

# Make a linear model based on data from samples with >200k reads where the rate of increase seems more-or-less linear (see plot2).
mod1<-lm(perc_compl~reads,data=rarefaction_dat[which(rarefaction_dat$reads>thresh),])
predframe<-data.frame("reads"=seq(thresh/100000,40)*10^5)
predframe$perc_compl<-predict(mod1,newdata=predframe)

plot1_plus1<-ggplot(data=rarefaction_dat)+
  geom_line(aes(x=reads/100000,y=perc_compl))+
  geom_line(aes(x=reads/100000,y=sd_upr_perc_compl),linetype="dashed")+
  geom_line(aes(x=reads/100000,y=sd_lwr_perc_compl),linetype="dashed")+
  geom_line(aes(x=reads/100000,y=max_perc_compl),linetype="dashed",colour="red")+
  geom_hline(yintercept = busco_score,linetype="dotdash")+ # This line is only relevant for my data, I'm
  geom_label(aes(x=Inf,y=busco_score,hjust=1,label="IsoSeq3 Clusters"))+ # This line is only relevant for my data, I'm
  geom_hline(yintercept = sp_g_busco,linetype="dotdash")+ # This line is only relevant for my data, I'm
  geom_label(aes(x=-Inf,y=sp_g_busco,hjust=0,label="Genome Assembly"))+ # This line is only relevant for my data, I'm
  ylim(0,100)+
  #  xlim(0,40)+
  #  geom_point(aes(x=size,y=mean,colour=data))+
  xlab(expression(paste("# Subsampled FL CCS Reads (x 10"^5,")")))+
  ylab("% Complete BUSCO Genes")+
  ggtitle(label=paste(dataset,"-",sampl))+
  my.theme+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))+
  geom_line(data=predframe[which(predframe$perc_compl<=100),],
            aes(x=reads/100000,y=perc_compl),colour="grey50",linetype="dashed")

ggsave(filename=here("figures",paste(dataset,"_plot1_plus1.png",sep="")),device = "png",dpi = 300,width = 15,height=10,units="cm")
plot1_plus1
dev.off()
ggsave(filename=here("figures",paste(dataset,"_plot1_plus1.pdf",sep="")),device = "pdf",dpi = 300,width = 15,height=10,units="cm")
plot1_plus1
dev.off()



# Load the final busco table
busco_table2<-read.table(here(paste("data/",dataset,sep=""),"sample1_2.polished_busco_full_table.tsv"),
                        header=FALSE,sep="\t",fill = TRUE)
colnames(busco_table2)<-c("busco_id","status","transcript","length","length2")
# get the total nr buscos from the busco_table
busco_total2<-length(unique(busco_table2$busco_id))
# Subset results to only complete (Single copy) or complete (duplicated) buscos
busco_table2<-busco_table2[which(busco_table2$status=="Complete" | busco_table2$status=="Duplicated"),]
busco_score2<-(length(unique(busco_table2$busco_id))/busco_total)*100

# Load the final reads table
trans_read_counts<-read.table(here(paste("data/",dataset,sep=""),
                                   "sample1_2.polished.cluster_report.csv"),
                              header=TRUE,sep=",")
# Need to rename the transcripts in the trans_read_counts table to match the busco output 
# (busco doesn't like the "/" character in names)
trans_read_counts$cluster_id<-gsub("/","_",trans_read_counts$cluster_id)

# How many reads per transcript
count<-tapply(trans_read_counts$read_id,INDEX = list(trans_read_counts$cluster_id),length)
transcript<-names(count)
counts_dat<-data.frame("cluster_id"=transcript,"count"=unname(count))

# How many transcripts with count >= 2
nrow(counts_dat[counts_dat$count >=2,])
n_reads2<-nrow(trans_read_counts)

plot1_plus2<-plot1<-ggplot(data=rarefaction_dat)+
  geom_line(aes(x=reads/100000,y=perc_compl))+
  geom_line(aes(x=reads/100000,y=sd_upr_perc_compl),linetype="dashed")+
  geom_line(aes(x=reads/100000,y=sd_lwr_perc_compl),linetype="dashed")+
  geom_line(aes(x=reads/100000,y=max_perc_compl),linetype="dashed",colour="red")+
  
  geom_hline(yintercept = busco_score,linetype="dotdash")+ # This line is only relevant for my data, I'm
  geom_label(aes(x=Inf,y=busco_score,hjust=1,label="IsoSeq3 Clusters"))+ # This line is only relevant for my data, I'm
  geom_hline(yintercept = sp_g_busco,linetype="dotdash")+ # This line is only relevant for my data, I'm
  geom_label(aes(x=-Inf,y=sp_g_busco,hjust=0,label="Genome Assembly"))+ # This line is only relevant for my data, I'm
  geom_point(aes(x=n_reads2/100000,y=busco_score2),colour="red",size=3)+
  ylim(0,100)+
  #  xlim(0,40)+
  #  geom_point(aes(x=size,y=mean,colour=data))+
  xlab(expression(paste("# Subsampled FL CCS Reads (x 10"^5,")")))+
  ylab("% Complete BUSCO Genes")+
  ggtitle(label=paste(dataset,"-",sampl))+
  my.theme+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))+
  geom_line(data=predframe[which(predframe$perc_compl<=100),],
            aes(x=reads/100000,y=perc_compl),colour="grey50",linetype="dashed")

ggsave(filename=here("figures",paste(dataset,"_plot1_plus2.png",sep="")),device = "png",dpi = 300,width = 15,height=10,units="cm")
plot1_plus2
dev.off()
ggsave(filename=here("figures",paste(dataset,"_plot1_plus2.pdf",sep="")),device = "pdf",dpi = 300,width = 15,height=10,units="cm")
plot1_plus2
dev.off()




# Plot sample1 + sample1_2
sampl<-"sample1_2" # One of "sample1", "sample1_2"

# Load the busco results table
sample1_2_rarefaction_dat<-read.table(here(paste("data/",dataset,sep=""),
                                 paste(sampl,"_subsampling_rarefaction_dat.csv",sep="")),
                            header=TRUE,sep=",")

sample1_2_busco_table<-read.table(here(paste("data/",dataset,sep=""),
                             paste(sampl,".polished_busco_full_table.tsv",sep="")),
                        header=FALSE,sep="\t",fill = TRUE)
colnames(sample1_2_busco_table)<-c("busco_id","status","transcript","length","length2")
# get the total nr buscos from the busco_table
sample1_2_busco_total<-length(unique(sample1_2_busco_table$busco_id))
# Subset results to only complete (Single copy) or complete (duplicated) buscos
sample1_2_busco_table<-sample1_2_busco_table[
  which(sample1_2_busco_table$status=="Complete" | sample1_2_busco_table$status=="Duplicated"),]
sample1_2_busco_score<-(length(unique(sample1_2_busco_table$busco_id))/sample1_2_busco_total)*100


plot1_plus3<-ggplot()+
  geom_line(data=rarefaction_dat,
            aes(x=reads/100000,y=perc_compl),colour="red")+
  geom_line(data=sample1_2_rarefaction_dat,
            aes(x=reads/100000,y=perc_compl),colour="blue")+
  
  geom_hline(yintercept = busco_score,linetype="dotdash",colour="red")+ # This line is only relevant for my data, I'm
  geom_label(aes(x=Inf,y=busco_score,hjust=1,label="IsoSeq3 Clusters (sample1)"))+ # This line is only relevant for my data, I'm
  
  geom_hline(yintercept = sample1_2_busco_score,linetype="dotdash",colour="blue")+ # This line is only relevant for my data, I'm
  geom_label(aes(x=Inf,y=sample1_2_busco_score,hjust=1,label="IsoSeq3 Clusters (sample1+2)"))+ # This line is only relevant for my data, I'm
  
  geom_hline(yintercept = sp_g_busco,linetype="dotdash")+ # This line is only relevant for my data, I'm
  geom_label(aes(x=-Inf,y=sp_g_busco,hjust=0,label="Genome Assembly"))+ # This line is only relevant for my data, I'm
  
  geom_point(aes(x=n_reads2/100000,y=busco_score2),colour="red",size=3)+

  geom_line(data=predframe[which(predframe$perc_compl<=100),],
            aes(x=reads/100000,y=perc_compl),colour="grey50",linetype="dashed")+

  ylim(0,100)+
  #  xlim(0,40)+
  #  geom_point(aes(x=size,y=mean,colour=data))+
  xlab(expression(paste("# Subsampled FL CCS Reads (x 10"^5,")")))+
  ylab("% Complete BUSCO Genes")+
  ggtitle(label=paste(dataset,"- sample1 + ",sampl))+
  my.theme+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))

ggsave(filename=here("figures",paste(dataset,"_plot1_plus3.png",sep="")),device = "png",dpi = 300,width = 15,height=10,units="cm")
plot1_plus3
dev.off()
ggsave(filename=here("figures",paste(dataset,"_plot1_plus3.pdf",sep="")),device = "pdf",dpi = 300,width = 15,height=10,units="cm")
plot1_plus3
dev.off()


