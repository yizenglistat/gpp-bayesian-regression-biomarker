#!/bin/sh
mkdir -p log
rm -rf log/*
git clone https://github.com/yizenglistat/gpp-bayesian-regression-biomarker.git
rm -rf output/ R/ run.r readme.md
cp -r gpp-bayesian-regression-biomarker/* .
rm -rf gpp-bayesian-regression-biomarker
# check it
# > ls
# update.sh R/ output/ run.r readme.md