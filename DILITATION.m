%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : DILITATION.m
%  Purpose : Computes the dilatation (areal strain) of a 2D displacement
%            field on scattered stations. The displacement components are
%            interpolated with RBFs (RBF.m) and the numerical dilatation
%            (1/2)*(du/dx + dv/dy) is compared with the analytical one.
%            Application: geodetic deformation analysis of a network.
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
close all
clc
format long g

m = 15 ; % m*m is the number of fixed points

x = linspace (-5, 5, m) ;
y = linspace (-5, 5, m) ;

[X, Y] = meshgrid (x, y) ;

X = X + randn (m, m) ;
Y = Y + randn (m, m) ;

U = sin (sqrt (X .^ 2 + Y .^ 2)) ;
V = cos (sqrt (X .^ 2 + Y .^ 2)) ;

x = reshape (X, [m^2 1]) ;
y = reshape (Y, [m^2 1]) ;
u = reshape (U, [m^2 1]) ;
v = reshape (V, [m^2 1]) ;

n = 60 ; % n*n is the number of the interpolated points

x0 = linspace (-4, 4, n) ;
y0 = linspace (-4, 4, n) ;

[X0, Y0] = meshgrid (x0, y0) ;

% U0 = sin (sqrt (X0 .^ 2 + Y0 .^ 2)) ;
% V0 = cos (sqrt (X0 .^ 2 + Y0 .^ 2)) ;

x0 = reshape (X0, [n^2 1]) ;
y0 = reshape (Y0, [n^2 1]) ;

du_dx =  (x0 ./ sqrt (x0 .^ 2 + y0 .^ 2)) .* cos (sqrt (x0 .^ 2 + y0 .^ 2)) ;
dv_dy = -(y0 ./ sqrt (x0 .^ 2 + y0 .^ 2)) .* sin (sqrt (x0 .^ 2 + y0 .^ 2)) ;

AD = (1/2) .* (du_dx + dv_dy) ;

h = waitbar (0, 'Please wait...') ; steps = n*n ;

tic
for i = 1:n*n
    [u0(i, :), du0_dx(i, :), du0_dy] = RBF (x, y, u, x0(i, :), y0(i, :), 10, 10, 'G') ;
    [v0(i, :), dv0_dx, dv0_dy(i, :)] = RBF (x, y, v, x0(i, :), y0(i, :), 10, 10, 'G') ;
    waitbar(i / steps)
end
toc

close (h)

ND = (1/2) .* (du0_dx + dv0_dy) ;

AD = reshape (AD, [n n]) ;
ND = reshape (ND, [n n]) ;

figure (1)
plot (x, y, 'r^'), hold on, quiver (x, y, u, v), xlabel ('X'), ylabel ('Y')
title ('Displacement'), legend ('Stations', 'Displacement'), grid on, axis equal

figure (2)
subplot (1,2,1)
imshow (AD), title ('Analytical Dilation')
subplot (1,2,2)
imshow (ND), title ('Numerical Dilation')

