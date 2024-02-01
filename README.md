Files and functions to compute square root of semi-infinite quasi-Toeplitz $M$-matrices
by Hongjia Chen, Hyun-Min Kim and Jie Meng
Feb. 1 2024


interpari.m
% function [gm,gp]=interpari(am1m,am1p)
% compute the coefficients of b(z), where b(z) satisfies b(z)^2-2b(z)+1=a(z)

sroot.m
% function [G,r] = sroot(B,X0,beta)
% compute B such that \beta(I-B)^2=A by fixed-point iteration

crsroot.m
% function [G,r] = crsroot(A)
% Compute the square root of quasi-Toeplitz M-matrix A by Cyclic reduction

 sroot_fixed1.m
% function [Cor_E,r] = sroot_fixed1( A, gn, gp )
% Compute the correction part of square root of M-matrix A=Ta+Cor_E, by fixed-point iteration


sda.m
% function [Cor_E,r] = sda(A,IN)
% Compute the correction part of square root of M-matrix A by structure-preserving algorithm

sroot_fixed1_finite.m
function [Cor_E,r] = sroot_fixed1_finite( A, Tb, N )
% Compute the correction part of the finite matrix A_N by fixed-point iteration
% A_N is the N by N principle section of A
%Tb is the Toeplitz part of the square root of A


sda_finite.m
% function [Cor_E,r] = sda_finite(A,INN,N)
% Compute the correction part of square root of A_N by structure-preserving algorithm
% A_N is the N by N principle section of A
% INN is the chosen initail point
