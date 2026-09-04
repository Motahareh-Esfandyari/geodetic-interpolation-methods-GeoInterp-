%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : chebyshev.m
%  Purpose : Function approximation using Chebyshev polynomials of the
%            first kind. Computes the Chebyshev coefficients and expands
%            the symbolic approximating polynomial.
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g

x0 = [0.1 0.6 1.1 1.6 2.1]' ; f = [1.1052 1.8221 3.0042 4.9530 8.1662]' ;

n = size (x0,1) ;

syms x

T (1 , :) = sym (1) ; T (2 , :) = x ;

for i = 2:(n - 1)
    
    T (i + 1 , :) = 2 * x * T (i , :) - T (i - 1 , :) ;
    
end

c (1) = sym((2 / n) * sum (f)) ;

for i = 2:n
    
    c (i , :) = sym ((2 / n) * sum (f .* subs (T (i) , x , x0))) ; 
    
end

p = (1 / 2) * c (1) + sum (c (2:end) .* T (2:end)) ;

p = expand (p) 