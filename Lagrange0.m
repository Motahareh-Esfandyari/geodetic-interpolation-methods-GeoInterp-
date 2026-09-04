%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Function: Lagrange0.m   [INCOMPLETE - WORK IN PROGRESS]
%  Purpose : Intended: Lagrange interpolation of satellite coordinates for
%            OrbitInterpLagrange.m. The implementation is unfinished.
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [PInt] = Lagrange0(xSat,ySat,zSat)
%LAGRANGE0 Summary of this function goes here
%   Detailed explanation goes here
for i = 1:size(xSat,1)

 l(i, j) = (x- x0(j))/(x0(i)-x0(j));
end