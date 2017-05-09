#!/usr/bin/Rscript
#Usage:
#Rscript buildRandomForest.r <dataset_training.csv>

library(randomForest)

# Reading arguments
args<-commandArgs(TRUE)
dataset_training_name = args[1]

nodesize_values = c(10, 20, 50, 100)
mtry_values = c(5, 10, 15)

# Reading dataset
dat_training = read.table(dataset_training_name, sep=";", header=T)
rsq <- list()
mse <- list()

for ( nodesize_v in nodesize_values ) {
	print(nodesize_v)
	rsq_i <- numeric()
	mse_i <- numeric()
	for ( mtry_v in mtry_values ) {
		print(mtry_v)
		rf = randomForest(MFE_frontier~., data=dat_training, nodesize=nodesize_v, ntree=500, mtry=mtry_v)
		rsq_i <- c(rsq_i, tail(rf['rsq'][[1]], 1))
		mse_i <- c(mse_i, tail(rf['mse'][[1]], 1))
		name = paste("RandomForest", nodesize_v, mtry_v, '.rsave', sep="-")
		save(rf, file=name)

	}
	rsq[[length(rsq)+1]] = rsq_i
	mse[[length(mse)+1]] = mse_i
}

print (rsq)


print (mse)
# Writing result
write.table(rsq,file="RandomForestResults-rsq.csv",sep=",",row.names=F)
#save(rsq, file="RandomForestResults-rsq.rsave")

write.table(mse,file="RandomForestResults-mse.csv",sep=",",row.names=F)
#save(mse, file="RandomForestResults-mse.rsave")

warnings()
