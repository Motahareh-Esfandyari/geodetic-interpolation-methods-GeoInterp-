%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : IDW.m
%  Purpose : Inverse Distance Weighting (IDW) interpolation experiment on
%            a noisy 2D grid. Studies how the interpolation error changes
%            with the number of neighbours (support domain size, k-NN).
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g

x = linspace(-5, 5, 20);
y = linspace(-5, 5, 20);

[X, Y] = meshgrid(x, y);

Z = (X-5).^2+(Y-5).^2;

X2 = X+randn(20, 20);
Y2 = Y+randn(20, 20);
Z2 = (X2-5).^2+(Y2-5).^2;

k = 1;
x = [];
y = []; 
z = [];

x = reshape(X2, [], 1) ;
y = reshape(Y2, [], 1) ;
z = reshape(Z2, [], 1) ;

x0 = x(255);
y0 = y(255);
z0 = z(255);

x(255,:) = [];
y(255,:) = [];
z(255,:) = [];

h = waitbar(0, 'Please wait...');

i = 1;

for n = 5:25
    
    [IDX,D] = knnsearch([x y],[x0 y0],'K',n) ;
    
    w = 1./D.^2;
    
    I(i,:) = sum(w' .*z(IDX))./sum(w');
        
    error = norm(I-z0);
    
    ERROR(i,:) = abs(I(i,:)-z0);
    
    Name = ['Error = ' , num2str(ERROR(i,:)) , ' & Number of Fixed points in neighbourhood = ' , num2str(n)];
    
    figure(i)
    subplot(2,2,1)
    plot(x,y,'ro') , hold on , plot(x(IDX),y(IDX),'bo') , hold on , plot(x0,y0,'m^')
    xlabel('X') , ylabel('Y') , title([Name , ' (2D)'])
    legend('All Points','Neighbourhood Points','Interpolated Point') , grid on
    
    subplot(2,2,2)
    plot(x(IDX),y(IDX),'bo') , hold on , plot(x0,y0,'m^')
    xlabel('X') , ylabel('Y') , title('Support Domain (2D)')
    legend('Neighbourhood Points','Interpolated Point') , grid on
    
    subplot (2,2,3)
    mesh(X,Y,Z) , hold on , plot3(x(IDX),y(IDX),z(IDX),'bo') , hold on
    plot3(x0,y0,z0,'g^') , hold on , plot3(x0,y0,I (i,:),'mo')
    xlabel('X') , ylabel('Y') , zlabel('Z') , title([Name , ' (3D)'])
    legend('All Points','Neighbourhood Points','Main Point','Interpolated Point') , grid on
    
    subplot(2,2,4)
    plot3(x(IDX),y(IDX),z(IDX),'bo') , hold on
    plot3(x0,y0,z0,'g^') , hold on , plot3(x0,y0,I (i,:),'mo')
    xlabel('X') , ylabel('Y') , zlabel('Z') , title('Support Domain (3D)')
    legend('Neighbourhood Points','Main Point','Interpolated Point') , grid on
    
    waitbar(i/21)
    
    i = i+1;
    
end

close(h)

figure(i)

n = 5:25;

plot(n, ERROR), xlabel('Number of fixed point in the support domain')
ylabel('Error of interpolation'),
title('Interpolation error according to the size of support domain')