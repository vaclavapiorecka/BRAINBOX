function [ gm ] = calculategm( x, m, not )
% The function is used to calculate spherical spline interpolation in 3D space.
%
% INPUTS:
%   x - vector of vlaues, -1 <= x <= 1
%   m - degree of Legendre polynomial
%   not - number of terms, degree of Legendre polynomial
%
% OUTPUTS:
%   gm - Legendre polynomial of x
%
% EXAMPLE:
%
% EXPLANATION:
%
% SEE ALSO:
%
% Author: Vaclava Piorecka (vaclava.piorecka@fbmi.cvut.cz, vaclava.piorecka@nudz.cz)
% Date:   
% 2017-10-01    creation of function

gm = 0;                     % alocation

for n = 1 : 1 : not         % create Legrende polynom
    Pn = legendre(n,x);
    gm = gm + ((2*n+1)/(n^m*(n+1)^m))*Pn(1) ;
end
gm = gm/(4*pi);

end

