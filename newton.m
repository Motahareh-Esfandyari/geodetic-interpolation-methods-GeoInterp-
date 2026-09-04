%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : newton.m
%  Purpose : Polynomial interpolation using Newton's divided differences.
%            Builds the symbolic interpolating polynomial for a given set
%            of sample points (x0, f).
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g

x0 = [0.1 0.6 1.1 1.6 2.1]' ;
f = [1.1052 1.8221 3.0042 4.9530 8.1662]' ;

n = size(x0,1);
m = n-1;
k = 1;

x2 = x0;
x1 = x0;
f0 = f;

syms x

for i = 1:m
        
	x1((n-(i-1)):end , :) = [];    
	x2(1:i , :) = [];
        
	N(i).Df = diff(f0)./(x2 - x1);
    
    dx = sym(1);
    
    for k = 1:i
        
        dx = dx*(x-x0(k));
        
    end
    
    N(i).Dx = dx;
        
	f0 = N(i).Df;
    x2 = x0;
    x1 = x0;
        
end

p = f(1);

for i = 1:m
    
    p = p+N(i).Df(1)*N(i).Dx;
    
end

p = expand(p) 
    
