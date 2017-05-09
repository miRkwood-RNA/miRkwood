#!/usr/bin/Rscript
#Usage:
#Rscript predictWithRandomForest.r <dataset_test.csv> <RandomForest.rsave>

library(randomForest)

# Reading arguments
args<-commandArgs(TRUE)
dataset_test_name = args[1]
rf_save_name = args[2]

# Loading objects from file
load(rf_save_name)
dat_test = read.table(dataset_test_name, sep=";", header=T)

# Prediction using RandomForest
p = predict(rf, dat_test)
dat_test['MFE_frontier_computed'] = p

# Writing result

rf_name_tmp = substr(rf_save_name, 13, nchar(rf_save_name) - 6)
output_name = paste('predicted', rf_name_tmp, '-', dataset_test_name, sep="")

write.csv(dat_test, file = output_name)
