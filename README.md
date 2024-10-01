## HEDE

This is the official repo for the paper [HEDE: Heritability estimation in high dimensions by Ensembling Debiased Estimators](https://arxiv.org/abs/2406.11184).

The first function in ```estimators.R``` is the HEDE function, which takes the data (X, y) and some optional parameters, and output a list of heritability estimates for different thresholding values. Many of these values will be the same for a certain dataset.

The remaining functions are other heritability estimators from relevant works.
