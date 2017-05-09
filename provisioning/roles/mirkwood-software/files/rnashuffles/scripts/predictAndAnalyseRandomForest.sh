#*/bin/sh
for NODESIZE in 2 5 8 10 15 20 40 50 60 80 100
do
	for MTRY in 7 10 13
	do
		RSAVE="RandomForest-$NODESIZE-$MTRY-.rsave"
		Rscript ~/shuffles/scripts/predictWithRandomForest.r dataset_test-50K.csv $RSAVE
		python ~/shuffles/scripts/evaluate_prediction_results.py --predict predicted-$NODESIZE-$MTRY--dataset_test-50K.csv --shuf dataset_test-50K.shuf > cmp-$NODESIZE-$MTRY.txt
	done
done
