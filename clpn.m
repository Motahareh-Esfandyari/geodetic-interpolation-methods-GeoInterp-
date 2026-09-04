%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Function: clpn.m  (formerly cp.m - renamed to match the function name)
%  Purpose : Computes Legendre polynomials Pn(z) and their derivatives
%            Pn'(z) for a complex argument z, by upward recurrence.
%            Used in geodesy for spherical-harmonic computations.
%  Note    : MATLAB conversion of the CLPN routine from Jin & Zhang,
%            "Computation of Special Functions" (Fortran original).
%  ------------------------------------------------------------------------
%  Adapted by : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [n,x,y,cpn,cpd]=clpn(n,x,y,cpn,cpd,varargin);

z=complex(x,y);
cpn(0+1)=complex(1.0d0,0.0d0);
cpn(1+1)=z;
cpd(0+1)=complex(0.0d0,0.0d0);
cpd(1+1)=complex(1.0d0,0.0d0);
cp0=complex(1.0d0,0.0d0);
cp1=z;
for  k=2:n;
cpf=(2.0d0.*k-1.0d0)./k.*z.*cp1-(k-1.0d0)./k.*cp0;
cpn(k+1)=cpf;
if(abs(x)== 1.0d0&y == 0.0d0);
cpd(k+1)=0.5d0.*x.^(k).*k.*(k+1.0d0);
else;
cpd(k+1)=k.*(cp1-z.*cpf)./(1.0d0-z.*z);
end;
cp0=cp1;
cp1=cpf;
end;  k=fix(n)+1;
return;
end

