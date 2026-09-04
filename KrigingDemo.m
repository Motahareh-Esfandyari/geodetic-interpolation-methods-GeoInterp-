%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : KrigingDemo.m  (formerly kriging0.m)
%  Purpose : Ordinary kriging demonstration. Generates a spatially
%            correlated random field, samples it at 500 random locations,
%            computes the experimental variogram (variogram.m), fits a
%            stable variogram model (variogramfit, from MATLAB FEX), and
%            interpolates the field by ordinary kriging (kriging.m) with
%            the kriging variance map.
%  Note    : Requires variogramfit.m from the MATLAB File Exchange.
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g


[X,Y] = meshgrid(0:500);
Z = randn(size(X));
Z = imfilter(Z,fspecial('gaussian',[40 40],8));

n = 500;
x = rand(n,1)*500;
y = rand(n,1)*500;
z = interp2(X,Y,Z,x,y);

subplot(2,2,1)
imagesc(X(1,:),Y(:,1),Z); axis image; axis xy
hold on
plot(x,y,'.k')
title('random field with sampling locations')

v = variogram([x y],z,'plotit',false,'maxdist',100);

subplot(2,2,2)
[dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable');
title('variogram')

[Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y);
subplot(2,2,3)
imagesc(X(1,:),Y(:,1),Zhat); axis image; axis xy
title('kriging predictions')
subplot(2,2,4)
contour(X,Y,Zvar); axis image
title('kriging variance')





