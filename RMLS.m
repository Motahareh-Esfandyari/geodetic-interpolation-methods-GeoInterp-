%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : RMLS.m
%  Purpose : Recursive (sequential) Moving Least Squares interpolation
%            with a full quadratic basis. Points are added to the support
%            domain one by one and the solution is updated recursively
%            until the coefficients converge; convergence of each
%            coefficient is plotted against the support domain size.
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
x = reshape(X2, [], 1);
y = reshape(Y2, [], 1);
z = reshape(Z2, [], 1);

x0 = x(100); 
y0 = y(100);
z0 = z(100);
x(100,:) = [];
y(100,:) = [];
z(100,:) = [];

[IDX,D] = knnsearch([x y],[x0 y0],'K',6);

A0 = [x(IDX).^2 y(IDX).^2 x(IDX).*y(IDX) x(IDX) y(IDX) ones(size(x(IDX),1),1)];
f0 = z(IDX);
Q0 = eye(size(f0,1));

xcap0 = inv(A0' *inv(Q0)*A0)*A0' *inv(Q0)*f0;
Qxcap0 = inv(A0' *inv(Q0)*A0);
error = 1;
k = 7;
i = 1;
h = waitbar(0, 'please wait');

while error > 1e-15
    
    [IDX,D] = knnsearch ([x y],[x0 y0],'K',k);
    
    idx = IDX(end);
    
    Ak = [x(idx).^ 2 y(idx).^2 x(idx).*y(idx) x(idx) y(idx) 1];
    fk = z(idx);
    Qk = eye(1);
    
    vk = fk-Ak*xcap0;
    Qxcapk = inv(inv(Qxcap0)+Ak' *inv(Qk)*Ak);
    xcapk = xcap0+Qxcapk*Ak' *inv(Qk)*vk;
    dx = xcapk-xcap0;
    error = norm(xcapk-xcap0); 
    Qxcap0 = Qxcapk; 
    xcap0 = xcapk;
    
    k = k+1;
    
    ERROR(i, :) = dx';
    i = i+1;
    w1 = norm(xcap0); 
    w2 = norm(xcapk);
    w = w2/w1;
    
    if w2 > w1
        w = w1/w2;
    end
    
    waitbar(w)
            
end

close (h)

a = xcapk(1);
b = xcapk(2);
c = xcapk(3);
d = xcapk(4);
e = xcapk(5);
f = xcapk(6);

zInt = a*x0^2+b*y0^2+c*x0*y0+d*x0+e*y0+f;

n = 7:(k-1);

figure(1)

subplot(2, 3, 1)
plot(n, ERROR(:, 1))
ylabel('Error of coefficient " a "')
xlabel('Number of fixed point in the support domain')
title('P(x,y) = ax^2 + by^2 + cxy + dx + ey + f')

subplot(2, 3, 2)
plot(n, ERROR(:, 2))
ylabel('Error of coefficient " b "')
xlabel('Number of fixed point in the support domain')
title('P(x,y) = ax^2 + by^2 + cxy + dx + ey + f')

subplot(2, 3, 3)
plot(n, ERROR(:, 3))
ylabel('Error of coefficient " c "')
xlabel('Number of fixed point in the support domain')
title('P(x,y) = ax^2 + by^2 + cxy + dx + ey + f')

subplot(2, 3, 4)
plot(n, ERROR(:, 4))
ylabel('Error of coefficient " d "')
xlabel('Number of fixed point in the support domain')
title('P(x,y) = ax^2 + by^2 + cxy + dx + ey + f')

subplot(2, 3, 5)
plot(n, ERROR(:, 5))
ylabel('Error of coefficient " e "')
xlabel('Number of fixed point in the support domain')
title('P(x,y) = ax^2 + by^2 + cxy + dx + ey + f')

subplot(2, 3, 6)
plot(n, ERROR(:, 6))
ylabel('Error of coefficient " f "')
xlabel('Number of fixed point in the support domain')
title('P(x,y) = ax^2 + by^2 + cxy + dx + ey + f')

