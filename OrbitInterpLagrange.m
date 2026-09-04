%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : OrbitInterpLagrange.m  (formerly FINAL.m - Lagrange variant)
%  Purpose : GPS satellite orbit interpolation from a precise ephemeris
%            (SP3) file using Lagrange interpolation (Lagrange0.m) instead
%            of cubic splines. Companion to OrbitInterpolation.m.
%  WARNING : Work in progress - Lagrange0.m is incomplete and the variable
%            names yInt/pInt are inconsistent in the plotting section.
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g

[sph,sp] = sp3Cread('igs13730.sp3');
sat = input('Enter the satellite number = ');
time = input('Enter the time [hour minute second] = ')*[3600 60 1]';

global OMEG
OMEG = 7.2921151467*10^(-5) ;
dp = sp ;
for i = 1:96
    
    t = (i-1)*900;
    R = [cos(OMEG*t) -sin(OMEG*t) 0 ; sin(OMEG*t)  cos(OMEG*t) 0 ; 0 0 1];
    dp(i).CRD = R*dp(i).CRD;
    
end

for i = 1:96
    
    epoch(i,:) = sp(i).Time;
    
    xSat(i,:) = dp(i).CRD(1,sat)*1000;
    ySat(i,:) = dp(i).CRD(2,sat)*1000;
    zSat(i,:) = dp(i).CRD(3,sat)*1000;
    
end

t = epoch(:,4).*3600+epoch(:,5).*60+epoch(:,6);

[xInt] = Lagrange0(t,xSat,time) ;
[pInt] = Lagrange0(t,ySat,time) ;
[zInt] = Lagrange0(t,zSat,time) ;

subplot(2,2,1)
plot3(xSat,ySat,zSat), title ('Satellite Orbit in an Inertial System')
xlabel('X[m]') , ylabel('Y[m]') , zlabel('Z[m]') , hold on
plot3(xInt,yInt,zInt,'ro')

subplot(2,2,2)
plot(t,xSat) , title('X component of orbit according to time')
xlabel('Time [second]') , ylabel('X [m]') , hold on
plot(time,xInt,'ro')

subplot(2,2,3)
plot(t,ySat) , title('Y component of orbit according to time')
xlabel('Time[second]') , ylabel('Y[m]') , hold on
plot (time,yInt,'ro')

subplot(2,2,4)
plot(t,zSat) , title('Z component of orbit according to time')
xlabel('Time [second]') , ylabel('Z [m]') , hold on
plot(time,zInt,'ro')