% read csv formatted matrix
M = csvread('sc-matrix.csv', 1, 1)

% test for sparsity (1 if sparse, 0 otherwise)
TF = issparse(M)