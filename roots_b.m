function r=roots_b(x)
% function r=myroots_b(a)
% r(i) is the solution of modulus less than 1 of equation b^2-2b+1-x=0;
% a,b,c column vectors 

a=1;b=-2;c=1-x;


d = sqrt(b.^2-4*a.*c);
a2 = 2*a;
x1 = (-b+d)./a2;
x2 = (-b-d)./a2;
if abs(x1)<1
    r=x1;
else
    r=x2;
end
