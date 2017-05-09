#!/bin/sh
#Usage: ./runRandomForest.sh <TRAINING.shuf> <TEST.shuf> <nodesize> <mtry> <ntree>
[ "$#" -eq 5 ] || die "Please provide parameters. Usage is ./runRandomForest.sh <TRAINING.shuf> <TEST.shuf> <nodesize> <mtry> <ntree>"
TRAINING=$1
TEST=$2
NODESIZE=$3
MTRY=$4
NTREE=$5
BASEDIR=$(dirname $0)
echo "Converting results file into CSV..."
python $BASEDIR/result2CSV.py --pvalue 1 $TRAINING
python $BASEDIR/result2CSV.py --pvalue 1 $TEST

#echo "Learning RandomForest..."
#Rscript $BASEDIR/buildRandomForest.r $TRAINING "RandomForest-30K.rsave"

echo "Predicting with RandomForest..."
Rscript $BASEDIR/predictWithRandomForest.r $TEST "RandomForest-30K.rsave"

echo "Evaluating prediction results quality..."
python $BASEDIR/evaluate_prediction_results.py --predict predicted_dataset.csv --shuf dataset_test.shuf

echo "Done."

