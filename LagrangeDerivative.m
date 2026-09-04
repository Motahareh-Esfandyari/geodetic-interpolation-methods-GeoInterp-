%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : LagrangeDerivative.m  (formerly darun_yabi_lagranj.m)
%  Purpose : Numerical interpolation/evaluation of a tabulated function
%            (here f = log(x)) at a point 'a' using Lagrange basis
%            coefficients computed from divided products of the nodes.
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

% nabejaee
clear
clc
close all

b = [0.4 0.5 0.7 0.8];      %CHANGE   % Matrix of Table
a = 0.6;                    %CHANGE   % Value

[m n]=size(b);
z=b;

for i=1:n
    x=z(1,i);
    s=log(x);               %CHANGE   % Function
    c(1,i)=s;   
end
for j=1:n
   for i=1:n
      z(i+1,j)=z(1,j)-z(1,i);
   end
end
y=z(2:n+1,:);
for j=1:n
    y(n+1,j)=1;
    for i=1:n
        if y(i,j)~=0
            y(n+1,j)=y(n+1,j)*y(i,j);
        end
    end
end
sum=1;
for i=1:n
    sum=sum*(a-z(1,i));
end
for i=1:n
    L(1,i)=((sum/(a-z(1,i)))/y(n+1,i));
end
L
x=0;
for i=1:n
    x=L(1,i)*c(1,i)+x;
end
out = x;
out