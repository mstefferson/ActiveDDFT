function [v, vFt] = decayexp2d( es, ls, n1, l1, n2, l2 )
% Calcualte distances
dx1    = l1/n1;
dx2    = l2/n2;
x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1]; 
x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
[x2m, x1m] = meshgrid( x2, x1 );
r = sqrt( x1m .^ 2 + x2m .^ 2 );
% soft shoulder potential
v = es .* exp( - r ./ ls );
% FT
vFt = fftshift( fftn( v ) );
