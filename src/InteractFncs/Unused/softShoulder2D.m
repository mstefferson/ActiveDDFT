function [v, vFt] = softShoulder2d( epsilon, a, R, Rs, n1, l1, n2, l2 )
% Calcualte distances
dx1    = l1/n1;
dx2    = l2/n2;
x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1]; 
x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
[x2m, x1m] = meshgrid( x2, x1 );
r2 = x1m .^ 2 + x2m .^ 2;
r8 = r2 .^ 4;
R8 = R ^ 8;
Rs8 = Rs ^ 8;
% soft shoulder potential
v = epsilon .* (  exp( - r8 ./ R8 ) + a .* exp( -r8 ./ Rs8 ) );
% FT
vFt = fftshift( fftn( v ) );
