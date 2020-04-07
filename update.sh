#!/bin/sh
mkdir -p log
rm -rf log/*
git clone https://github.com/yizenglistat/gpp-bayesian-regression-group-testing.git
rm -rf output/ R/ run.r readme.md
cp -r gpp-bayesian-regression-group-testing/* .
rm -rf gpp-bayesian-regression-group-testing
rm -rf output/IT/*.csv output/MPT/*.csv output/AT/*.csv output/DT/*.csv
# check it
# > ls
# update.sh R/ output/ run.r readme.md