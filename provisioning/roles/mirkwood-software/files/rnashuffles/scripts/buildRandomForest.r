#!/usr/bin/Rscript
#Usage:
#Rscript buildRandomForest.r <nodesize> <mtry> <ntree> <dataset_training.csv> <RandomForest.rsave>

library(randomForest)

# Reading arguments
args<-commandArgs(TRUE)

nodesize_value = args[1]
mtry_value = args[2]
ntree_value = args[3]
dataset_training_name = args[4]
rf_save_name = args[5]

# Reading dataset
dat_training = read.table(dataset_training_name, sep=";", header=T)

# Learning with RandomForest
rf = randomForest(MFE_frontier~., data=dat_training, nodesize=nodesize_value, ntree=ntree_value, mtry=mtry_value)

# Writing result
save(rf, file=rf_save_name)
