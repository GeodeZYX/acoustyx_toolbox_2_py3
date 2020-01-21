#!/bin/bash

hosnam=`hostname`
#git add ./scripts/benchmark_rt/*
git commit -a -m $hosnam
git push origin master
