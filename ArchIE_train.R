#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
train <- data.frame(fread(args[1]))
LABEL_COL=209

lr <- train[,1:LABEL_COL] # cut down the extra columns -- 209 is the label (0=not archaic, 1=archaic)
model <- glm(V209 ~ .,family=binomial(link='logit'),data=lr) # train the model
save.image(args[2])