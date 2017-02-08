%  [mayerFnc, mayerFncFt] = mayerFncHs(n1, n2, n3, l1, l2, l3, rSpheres)
%
% Calculates mayer function and ft for spheres.
%
function [mayerFnc, mayerFncFt] = mayerFncHs(n1, n2, n3, l1, l2, l3, rSpheres)
% Allocate
mayerFnc = zeros(n1,n2,n3);
% Calcualte distances
dx1    = l1/n1;
dx2    = l2/n2;
dx3    = l3/n3;
epsilon = 0.00001;
tooFar =  4 .* rSpheres .* rSpheres + epsilon;
x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1]; 
x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
x3 = dx3 .* [ 0 : 1 : n3 / 2  -n3/2+1 : 1 : -1];
[x2m, x1m, x3m] = meshgrid( x2, x1, x3 );
r2 = x1m .^ 2 + x2m .^ 2 + x3m .^ 2;
% Mayer fnc
ind2Accept = r2 < tooFar;
mayerFnc( ind2Accept) = -1;
% FT
mayerFncFt = fftshift( fftn( mayerFnc ) );
end
