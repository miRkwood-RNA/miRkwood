#!/bin/sh
# Helper script to run the browser tests using Cucumber (though bundle)
# Argument is the target environment to test:
# LOCAL for localhost, VAGRANT for the Vagrant virtual machine
# BIOINFOTEST and BIOINFO for the Bonsai qualification and production servers.

TARGET=$1

if [ "$TARGET" = "LOCAL" ]; then
    ROOT_URL="http://127.0.0.1"
elif [ "$TARGET" = "VAGRANT" ]; then
    ROOT_URL="http://192.168.33.20"
elif [ "$TARGET" = "BIOINFOTEST" ]; then
    ROOT_URL="http://bioinfotest.lifl.fr"
elif [ "$TARGET" = "BIOINFO" ]; then
    ROOT_URL="http://bioinfo.lifl.fr"
else
    echo "No argument supplied for Target environment"
    echo "Usage: `basename $0` [ LOCAL | VAGRANT | BIOINFOTEST | BIOINFO ]"
    exit 1
fi

if [ "$CI_MODE" ]; then
    ADDITIONAL='--format pretty --format junit --out . --format html --out browsertests-report.html HEADLESS=True'
else
    ADDITIONAL='--format pretty'
fi

if [ "$NO_EXECUTION" ]; then
    ADDITIONAL="$ADDITIONAL --tags ~@execution"
fi

echo "Running browser tests from $ROOT_URL"
MIRKWOOD_URL="$ROOT_URL/cgi-bin/mirkwood/web_scripts/"
MIRKWOOD_HOME_URL="$ROOT_URL/mirkwood/"

bundle exec cucumber $ADDITIONAL MIRKWOOD_URL=$MIRKWOOD_URL MIRKWOOD_HOME_URL=$MIRKWOOD_HOME_URL
