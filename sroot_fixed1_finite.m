function [Cor_E,r] = sroot_fixed1_finite( A, Tb, N )
% Compute the correction part of the finite matrix A_N
% A_N is the N by N principle section of A
% A, Tb are  semi-infinite matrices
verb = true; maxit = 1000; epsi = 1.e-13; cqtoption( 'threshold', 10^(-15) );
D = ones( 1, N );
I = diag( D );

FTb = Tb( 1 : N, 1 : N );
Tbb = Tb^2;
T12T21 = Tbb( 1 : N, 1 : N ) - FTb^2;

A1 = I - A( 1 : N, 1 : N );

Q11 = A1 + T12T21;

r = zeros( maxit, 1 );
Q2 = 2 * I - FTb;
Q1 = I - FTb;
 X = zeros( N, N );

err = 1;
for k = 1 : maxit
    Xold = X;
    errold = err;
    X = ( Q2 - Xold )^(-1) * Q11 - FTb;
    err = norm( ( Q1 - X )^2 - ( I - Q11 ), inf ) / norm( I - Q11, inf );
  
    
     if verb
        fprintf( 'step=%d, err=%d\n', k, err ); 
     end
     r(k) = err;
    
if err < epsi || ( err >= errold && k > 1 ),break;end
 
end

Cor_E = X;

if ( k == maxit )
fprintf( 'Warning: reached the max number of iterations' );
end
