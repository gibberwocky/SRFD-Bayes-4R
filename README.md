# SRFD-Bayes-4R

This repository is a *rough* R port of the MATLAB code for the SRFD-Bayes method (https://github.com/Astaxanthin/SRFD-Bayes) described by Zhou et al. This code has not yet been optimized, so there is no doubt room for improvement - which I will endeavour work on when time allows.

Zhou X, Cheng Z, Dong M, *et al.* Tumor fractions deciphered from circulating cell-free DNA methylation for cancer early diagnosis. *Nature Communications*, 2022, **13**(1): 1-13.

## R
The R code has been written in R 4.2.3 and has the following dependencies:
* matlab
* R.matlab
* pracma
* MASS
* caret
* e1071
* klaR
* RSNSS

## Usage
The MATLAB workflow outlined in main.m from Astaxanthin's repository is broadly replicated in `SRFD-Bayes.R` while parameters are defined in `R/Initialization.R`.
