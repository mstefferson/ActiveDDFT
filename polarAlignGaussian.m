% polar aligned potential with no gaussian drop-off in position
function [v, vFt] = polarAlignGaussian( es1,  ls1, n1, n2, n3, l1, l2, l3 );
% calculate distances
dx1    = l1/n1;
dx2    = l2/n2;
x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1]; 
x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
[x2m, x1m] = meshgrid( x2, x1 );
r2 = x1m .^ 2 + x2m .^ 2;
% calulate angle
phi = l3/n3 * (1:n3-1);
% polar alignedment with gaussian dropoff
v = es1 .*  cos( phi ) .* exp( - r2 / (2 * ls1 ^2 ) );
% Ft
vFt = fftshift( fftn( v ) );
end


