%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Function: CubicSpline.m
%  Purpose : Natural cubic spline interpolation. Solves the tridiagonal
%            system for the second derivatives and evaluates the spline
%            at the query point x0.
%  Inputs  : x, y - sample points (sorted by x) ; x0 - query point
%  Outputs : yInt - interpolated value at x0
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [yInt] = CubicSpline(x,y,x0)


h = diff(x);
dy = diff(y);
b = dy./h;

for i = 2:size(h,1)
    
    U(i-1 , :) = 2.*(h(i - 1 , :)+h(i , :));
    V(i-1 , :) = 6.*(b(i , :)-b(i-1 , :));
    
end

n = size(U,1);
B = zeros(n,n); 
B(1,1:2) = [U(1) h(2)];
B(n,n - 1:n) = [h(n) U(n)];

for i = 2:n-1
    
    B(i , i - 1:i + 1) = [h(i) U(i) h(i + 1)];
    
end

z = inv(B)*V;
z = [0 ; z ; 0];
[r1,~] = find(x==x0); 

if isempty(r1)~=1
    
    yInt = y(r1,:);
    
else
    
    [r2,~] = find(x < x0);
    i = r2(end);

    yInt = z(i+1)/(6*h(i))*(x0-x(i))^3+z(i)/(6*h(i))* ...
        (x(i+1)-x0)^3+(y(i+1)/h(i)-(h(i)*z(i+1))/6)* ...
        (x0-x(i))+(y(i)/h(i)-(h(i)*z(i))/6)*(x(i+1)-x0);
    
end

end