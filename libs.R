# Clean workspace
rm(list=ls())
len = length

# Load libraries and user functions
# library(doSNOW)
# cl = makeCluster(2, type = 'SOCK')
# registerDoSNOW(cl)
library(Rcpp)
source('r_functions.R')
sourceCpp('cpp_functions.cpp')
library(gtools)
library(pracma)
library(extraDistr)
library(foreach)
library(microbenchmark)
library(progress)
library(seqinr)