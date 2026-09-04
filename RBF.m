%% ========================================================================
%  GeoInterp : Numerical Interpolation & Approximation Methods in Geodesy
%  ------------------------------------------------------------------------
%  Function: RBF.m
%  Purpose : Local Radial Basis Function (RBF) interpolation with
%            analytical first derivatives. Supports Gaussian ('G'),
%            Multiquadric ('MQ') and Inverse Multiquadric ('IMQ') kernels;
%            support domain selected by k-nearest neighbours.
%  Inputs  : x, y, u - scattered data ; x0, y0 - query point ;
%            n, m - support domain sizes ; BaseFunction - 'G'|'MQ'|'IMQ'
%  Outputs : u0 - interpolated value ; du0_dx, du0_dy - derivatives
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [u0, du0_dx, du0_dy] = RBF (x, y, u, x0, y0, n, m, BaseFunction)


[IDX, r] = knnsearch ([x, y], [x0, y0], 'k', n);

IDX = IDX' ;
r = r' ;
D = max (r);
N = n;
c = 1.25 .* D ./ sqrt (N);
c = c' ;

for i = 1:size (IDX, 2)
    xS (:, i) = x (IDX(:, i), :) ;
    yS (:, i) = y (IDX(:, i), :) ;
    uS (:, i) = u (IDX(:, i), :) ;
end

idx = [] ;
for i = 1:m
    for j = 1:n
        
        idx = [idx; [i, j]] ;
        
    end
end

rij = sqrt ((xS(idx (:, 2)) - xS(idx (:, 1))) .^ 2 +...
            (yS(idx (:, 2)) - yS(idx (:, 1))) .^ 2) ;
if strcmp (BaseFunction, 'G') == 1
    
    PHI = exp (-c .* rij .^ 2 ) ;
    PHI = reshape (PHI, [n, m]) ;
    
    L = (PHI' * PHI) ^ (-1) * PHI' * uS ;
    
    phi = exp (-c .* r (1:m, :) .^ 2) ;
    
    dp_dx = -2 .* c .* (x0 - xS (1:m, :)) .* exp (-c .* r (1:m, :) .^ 2) ;
    dp_dy = -2 .* c .* (y0 - yS (1:m, :)) .* exp (-c .* r (1:m, :) .^ 2) ;
    
elseif strcmp (BaseFunction, 'MQ') == 1
    
    PHI = sqrt (1 + (c .* rij) .^ 2) ;
    PHI = reshape (PHI, [n, m]) ;
    
    L = (PHI' * PHI) ^ (-1) * PHI' * uS ;
    
    phi = sqrt (1 + (c .* r (1:m, :)) .^ 2) ;
    
    dp_dx = (c .^ 2 .* (2 .* x0 - 2 .* xS (1:m, :))) ./ (2 .* (r (1:m, :) .^ 2 .* c .^ 2 + 1) .^ (1/2)) ;
    dp_dy = (c .^ 2 .* (2 .* y0 - 2 .* yS (1:m, :))) ./ (2 .* (r (1:m, :) .^ 2 .* c .^ 2 + 1) .^ (1/2)) ;
    
elseif strcmp (BaseFunction, 'IMQ') == 1
    
    PHI = 1 ./ sqrt (1 + (c .* rij) .^ 2) ;
    PHI = reshape (PHI, [n, m]) ;
    
    L = (PHI' * PHI) ^ (-1) * PHI' * uS ;
    
    phi = 1 ./ sqrt (1 + (c .* r (1:m, :)) .^ 2) ;
    
    dp_dx = -(c .^ 2 .* (2 .* x0 - 2 .* xS (1:m, :))) ./ (2 .* (r (1:m, :) .^2 .* c .^ 2 + 1) .^ (3/2)) ;
    dp_dy = -(c .^ 2 .* (2 .* y0 - 2 .* yS (1:m, :))) ./ (2 .* (r (1:m, :) .^2 .* c .^ 2 + 1) .^ (3/2)) ;
    
end

u0 = L' * phi ;

du0_dx = L' * dp_dx ;
du0_dy = L' * dp_dy ;

end