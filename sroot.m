function [G,r] = sroot(B,X0,beta)
% function G = sroot(A)
% A=beta * (I - B) with B\geq 0 and \|B\|_{\infty}<1
% Compute the square root of quasi-Toepltitz M-matrix A by Binomial
% iteration
% A^(1/2) = beta^(1/2) * (I - G)
verb = true; maxit = 1000; epsi = 1.e-13;  cqtoption('threshold',10^(-15))
I = cqt( 1, 1 );
X = X0;
err = 1;
r = zeros( maxit, 1 );
for k = 1:maxit
    Xold   = X;
    X      = 1/2 * (B + Xold^2);
    
    errold = err;
    err    = norm( (I-X)^2-(I-B), inf )/norm( I-B, inf );
   
    
    
     if verb
         fprintf( 'step=%d, err=%d\n', k, err ); 
    end
    
  r(k)   = err;
 
    
if err < epsi || (err - errold >= 0 && k > 1), break; end

end
G = beta^(1/2) * (I - X);
if (k == maxit)
fprintf( 'Warning: reached the max number of iterations' );
end
