#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggforce)
library(ggpubr)
library(ggplot2)

# dat_mnm <- read.delim("output/sprime/GutenkunstThreePop/topSegment_wMNM.out.score", header=F,sep='')
# dat_womnm <- read.delim("output/sprime/GutenkunstThreePop/topSegment_woMNM.out.score", header=F,sep='')
dat_mnm <- read.delim(args[1], header=F, sep='')
dat_womnm <- read.delim(args[2], header=F, sep='')

dat_mnm$type = "wMNM"
dat_womnm$type = "woMNM"
dat <- rbind(dat_mnm, dat_womnm)

# pdf("output/sprime/GutenkunstThreePop/sinaplot_topSegment_wMNM_vs_woMNM.out.score.pdf", width=8, height=6)
pdf(args[3], width=8, height=6)
ggplot(dat, aes(type, V8)) +
	geom_sina(aes(color = type), size = 1) +
	theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) +
	xlab("") + ylab("Top Sprime score per segment") +
	scale_color_manual(values =  c("#FC4E07", "#00AFBB")) +
	stat_compare_means(method="wilcox.test", comparisons = list(c("wMNM", "woMNM"))) +
	stat_compare_means(label.y = 50, label.x=1.5)
dev.off()

