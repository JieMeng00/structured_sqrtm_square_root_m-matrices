function [G,r] = crsroot(A)
% function G = crsroot(A); 
% Compute the square root of quasi-Toeplitz M-matrix A by Cyclic reduction
verb = true; maxit = 1000; epsi = 1.e-13;  cqtoption('threshold',10^(-15))
I = cqt( 1, 1 );
Y = I - A;
Z = 2 * ( I + A );
err= 1;
r = zeros( maxit, 1 );
for k=1:maxit
    Yold = Y;
    Zold = Z;
    
    Y    = -Yold * Zold^(-1) * Yold;
    Z    =  Zold + 2 * Y;
    
    errold = err;
    err    = norm( 1/16*Z^2 - A, inf )/norm( A, inf );
   
    
    
    if verb 
      fprintf( 'step=%d, err=%d\n', k, err ); 
    end
 
      r(k)   = err;   
      
    if err < epsi || ( err >= errold && k > 1 ), break; end
    
      
 
end
G = 1/4 * Z;

if ( k == maxit )
fprintf( 'Warning: reached the max number of iterations' );
end
