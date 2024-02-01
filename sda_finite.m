function [Cor_E,r] = sda_finite(A,INN,N)
% Compute the correction part of square root of A_N
% A_N is the N by N principle section of A
verb = true; maxit = 1000; epsi = 1.e-13; cqtoption( 'threshold', 10^(-15) );
D = ones( 1, N );
I = diag( D );
B = cqt( 1, 1 ) - A;
r = zeros( maxit, 1 );
Rbb = ( INN - 2 * cqt( 1, 1 ) ) * INN + B;

IN = INN( 1 : N, 1 : N );

Rb = Rbb( 1 : N, 1 : N );

Tbb = INN^2;
T12T21 = Tbb( 1 : N, 1 : N ) - IN^2;

FA = B( 1 : N, 1 : N ) + T12T21;

S = ( 2 * I - IN )^( -1 );
Q = S;
P = S * Rb;
E = IN + P;
F = S;
T = I - IN;
err = norm( ( T - P )^2 - ( I - FA ), inf ) / norm( ( I - FA ), inf );

for k = 1 : maxit       
    Pold = P;
    Qold = Q;
    Eold = E;
    Fold = F; 

    FPQ = Fold * ( I - Pold * Qold )^(-1);
      P = Pold + FPQ * Pold * Eold;

 errold = err;
    err = norm( ( T - P )^2 - ( I - FA ), inf ) / norm( I - FA, inf );
   
     if verb
         fprintf( 'step=%d, err=%d\n', k, err ); 
     end


   r(k) = err;  
    
if err<epsi || ( err >= errold && k > 1 ),break;end  

    EQP = Eold * ( I - Qold * Pold )^( -1 ); 
      F = FPQ * Fold;
      E = EQP * Eold;
      Q = Qold + EQP * Qold * Fold;
end

Cor_E = P; % the non-Toeplitz part

if ( k == maxit )
fprintf( 'Warning: reached the max number of iterations' );
end