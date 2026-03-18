#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(stringr))

options(scipen = 999)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

df = fread(args[1], header=TRUE, sep=',')
df$Sample <- str_remove(df$Sample_Name, "_toulligqc")

n <- nrow(df)
width <- 14 + log10(n) * 2

#N50
p1 <- ggplot(df, aes(x=Sample, y=N50)) +
        geom_bar(stat="identity", color="black", fill='steelblue') +
        theme_bw() +
        xlab("Sample") +
        ylab("N50") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16)) +
        labs(title="Per Sample N50 value")

ggsave(p1, filename=paste0(args[2], "N50.png"), height=12, width = width, units = "in", dpi=300)

#Yield
p2 <- ggplot(df, aes(x=Sample, y=Yield)) +
        geom_bar(stat="identity", color="black", fill='green') +
        theme_bw() +
        xlab("Sample") +
        ylab("Yield[G]") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16)) +
        labs(title="Per Sample Yield")

ggsave(p2, filename=paste0(args[2], "yield.png"), height=12, width = width, units = "in", dpi=300)

#Mean Phred
p3 <- ggplot(df, aes(x=Sample, y=Pass_Reads_QScore_Mean)) +
  geom_bar(stat="identity", color="black", fill='red') +
        theme_bw() +
        xlab("Sample") +
        ylab("Mean Phred Score") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16)) +
        labs(title="Per Sample Mean Phred Score") +
        geom_hline(yintercept=mean(df$Pass_Reads_QScore_Mean), linetype='dashed') +
        labs(caption='Dashed line indicates mean over all samples.')

ggsave(p3, filename=paste0(args[2], "phred.png"), height=12, width = width, units = "in", dpi=300)

#Passed reads in percent
p4 <- ggplot(df, aes(x=Sample, y=Read_Pass_Percent)) +
        geom_bar(stat='identity', color="black", fill="orange") +
        theme_bw() +
        xlab("Sample") +
        ylab("Passed reads[%]") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16)) +
        labs(title="Passed reads per sample")

ggsave(p4, filename=paste0(args[2], "passed_reads.png"), height=12, width = width, units = "in", dpi=300)

#Read Count
p5 <- ggplot(df, aes(x=Sample, y=Read_Count)) +
        geom_bar(stat='identity', color="black", fill="purple") +
        theme_bw() +
        xlab("Sample") +
        ylab("Read count") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16)) +
        labs(title="Read count per sample")

ggsave(p5, filename=paste0(args[2], "read_count.png"), height=12, width = width, units = "in", dpi=300)

#L50
p6 <- ggplot(df, aes(x=Sample, y=L50)) +
        geom_bar(stat='identity', color="black", fill="steelblue") +
        theme_bw() +
        xlab("Sample") +
        ylab("L50") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16)) +
        labs(title="Per Sample L50 value")

ggsave(p6, filename=paste0(args[2], "L50.png"), height=12, width = width, units = "in", dpi=300)

#Read Length
df2 <- df[,c(13,8,9,10)]
colnames(df2) <- c("Sample", "Mean_Length", "Min", "Max")
to.plot <- reshape2::melt(df2)
p7 <- ggplot(to.plot, aes(x=Sample, y=value, group=variable, color=variable)) +
        geom_point(size=2.5) +
        geom_line() +
        theme_bw() +
        xlab("Sample") +
        ylab("Read length [bp]") +
        labs(title="Mean read length per sample") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16)) +
        scale_y_log10()

ggsave(p7, filename=paste0(args[2], "read_length.png"), height=12, width = width, units = "in", dpi=300)

if(length(args) > 2){
  df2 <- data.frame()
  con <- file(args[3], "r")
  
  while(TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    
    print(paste0("Processing: ", line))
    
    sample <- fread(line, header=F, sep="\t")
    colnames(sample) <- c("chr","start","end","modbase","score","strand","start1","end1","color","Nvalid_cov","fraction_mod","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall")
  
    mean_cov <- mean(sample$Nvalid_cov)
    mean_meth <- mean((sample$fraction_mod / 100))
    
    sample_name <- basename(line)
    df2 <- rbind(df2, data.frame("Samplename"=sample_name, "Mean_Cov"=mean_cov, "Mean_Meth"=mean_meth))
  }
  
  df2$Sample <- str_remove(df2$Samplename, ".bed")
  
  #Mean Cov
  p8 <- ggplot(df2, aes(x=Sample, y=Mean_Cov)) +
    geom_bar(stat="identity", fill="brown") +
    theme_bw() +
    xlab("Sample") +
    ylab("Mean Coverage") +
    labs(title="Mean coverage per sample") +
    geom_hline(yintercept=mean(df2$Mean_Cov), linetype='dashed') +
    labs(caption='Dashed line indicates mean over all samples.') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16))
  
  ggsave(p8, filename=paste0(args[2], "mean_cov.png"), height=12, width = width, units = "in", dpi=300)
  
  #Mean Meth
  p9 <- ggplot(df2, aes(x=Sample, y=Mean_Meth)) +
    geom_bar(stat="identity", fill="violet") +
    theme_bw() +
    xlab("Sample") +
    ylab("Mean Methylation[0-1]") +
    labs(title="Mean methylation per sample") +
    ylim(c(0,1)) + 
    geom_hline(yintercept=mean(df2$Mean_Meth), linetype='dashed') +
    labs(caption='Dashed line indicates mean over all samples.') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), plot.title=element_text(hjust=0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16,vjust=-0.5),text=element_text(size=16))
  
  ggsave(p9, filename=paste0(args[2], "mean_meth.png"), height=12, width = width, units = "in", dpi=300)
  
  close(con)
}