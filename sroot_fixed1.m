function [Cor_E,r] = sroot_fixed1( A, gn, gp )
% function G = sroot(A)-Tb; Tb=cqt(gn,gp)
% Compute the correction part of square root of M-matrix A
verb = true; maxit = 1000; epsi = 1.0e-13; cqtoption( 'threshold', 10^(-15) );
 I = cqt( 1, 1 );
A1 = I - A;
Tb = cqt( gn, gp );
 r = zeros( maxit, 1 );
Q2 = 2 * I - Tb;
Q1 = I - Tb;
 X = cqt( 0, 0 );
err = 1;
for k = 1 : maxit
      Xold = X;
    errold = err;
         X = ( Q2 - Xold )^( -1 ) * A1 - Tb;
       err = norm( ( Q1 - X )^2 - A, inf ) / norm( A, inf );
 
    
     if verb
         fprintf( 'step=%d, err=%d\n', k, err ); 
     end

  r(k)   = err;  
    
if err < epsi || ( err >= errold && k > 1 ),break;end
 
end

Cor_E = X;

if ( k == maxit )
fprintf('Warning: reached the max number of iterations');
end
