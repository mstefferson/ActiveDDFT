%  [mayerFnc, mayerFncFt] = mayerFncHd(n1, n2, n3, l1, l2, rDisks)
%
% Calculates mayer function and ft for disks. If there is a third
% dimension, it will rep the mayer fnc in that dimension
% dDisks: diameter of disks
function [mayerFnc, mayerFncFt] = mayerFncHd(n1, n2, n3, l1, l2, dDisks)
% Allocate
mayerFnc = zeros(n1,n2);
% Calcualte distances
dx1    = l1/n1;
dx2    = l2/n2;
epsilon = 0.00001;
tooFar =  dDisks .* dDisks + epsilon;
x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1]; 
x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
[x2m, x1m] = meshgrid( x2, x1 );
r2 = x1m .^ 2 + x2m .^ 2;
% Mayer fnc
ind2Accept = r2 < tooFar;
mayerFnc( ind2Accept) = -1;
mayerFnc = repmat( mayerFnc, [1, 1, n3]);
% FT
mayerFncFt = fftshift( fftn( mayerFnc ) );
end
