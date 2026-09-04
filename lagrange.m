%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : lagrange.m
%  Purpose : Polynomial interpolation using Lagrange basis polynomials.
%            Builds the symbolic interpolating polynomial for a given set
%            of sample points (x0, f).
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
close all
clc
format long g

x0 = [1;4;6];
f = [0;1.386294;1.79176];
n = size(x0,1); 

syms x

l = sym(zeros(n, n));

for i = 1:n
    for j = 1:n
        if i == j
            l(i, j) = sym(1);
        else
            l(i, j) = (x- x0(j))/(x0(i)-x0(j));
        end
    end
end

for i = 1:n
    L(:, i) = prod(l(i, :));
end

P = expand(L*f) 