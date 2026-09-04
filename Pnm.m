%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Script  : Pnm.m
%  Purpose : Driver for clpn.m - evaluates the complex Legendre
%            polynomials Pn(z) and derivatives at z = 3 + 2i, n = 0..5,
%            and prints a formatted table.
%  WARNING : Lines after the header block contain leftovers from the
%            Fortran conversion that are not valid MATLAB and must be
%            removed or commented before running.
%  ------------------------------------------------------------------------
%  Adapted by : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g


z = 3.0 +2.0 i
 n = Re[Pn(z)]Im[Pn(z)]Re[Pn'(z)]Im[Pn'(z)]
   
n=[];x=[];
y=[];
cpn=[];
cpd=[];
cpn=zeros(1,100+1);
cpd=zeros(1,100+1);

fprintf(1,'%s \n','  please enter nmax, x and y(z=x+iy)');
%        READ(*,*)N,X,Y
n=5;
x=3.0;
y=2.0;
fprintf(1,[repmat(' ',1,3),'x =','%5.1g',',  ','y =','%5.1g' ' \n'],x,y);
fprintf(1,'%0.15g \n');
[n,x,y,cpn,cpd]=clpn(n,x,y,cpn,cpd);
fprintf(1,'%s ','  n    re[pn(z)]im[pn(z)]re[pn''(z)]');fprintf(1,'%s \n', '   im[pn''(z)]');
fprintf(1,'%s ',' ---------------------------------------------');fprintf(1,'%s \n', '--------------');
for  k=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%14.6g',1,4) ' \n'],k,cpn(k+1),cpd(k+1));
end;  k=n+1;

