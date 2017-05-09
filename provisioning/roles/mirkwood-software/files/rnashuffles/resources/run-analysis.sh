#!/bin/bash
BASEDIR=$(dirname $0)
pip install -r $BASEDIR/analysis-requirements.txt
pylint --rcfile=.pylintrc shuffles > pylint.txt
pep8 > pep8.txt
