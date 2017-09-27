% [gridObj] = GridMakerPBCxk(n1,n2,n3,l1,l2,l3)
%
% Makes the gridObj. Real and k space coordinate vectors
%
function [gridObj] = GridMakerPBCxk(n1,n2,n3,l1,l2,l3)
% Make pos grid spacings
dx1 = l1/n1;
dx2 = l2/n2;
dx3 = l3/n3;
% Make vectors and grids
gridObj.x1 = dx1 * ( -n1/2: n1/2 - 1 );
gridObj.x2 = dx2 * ( -n2/2: n2/2 - 1 );
gridObj.x3 = dx3 * ( -n3/2: n3/2 - 1 );
% Make k-space spacings
dk1 = 2*pi/l1;
dk2 = 2*pi/l2;
dk3 = 2*pi/l3; %  = 1
% Make k vectors and grids
k1_max = pi / dx1; %Maximum spatial k-vector allowed by grid
k2_max = pi / dx2; %Maximum spatial k-vector allowed by grid
k3_max = pi / dx3; %Maximum angular k-vector
% For generalization purposes, handle N = 1 cases differently.
if n1 == 1
  gridObj.k1 = 0;
else
  gridObj.k1 = ( -k1_max: dk1: (k1_max - dk1) );
end
if n2 == 1
  gridObj.k2 = 0;
else
  gridObj.k2 =  ( -k2_max: dk2: (k2_max - dk2) );
end
if n3 == 1
  gridObj.k3 = 0;
else
  gridObj.k3 = ( -k3_max: dk3: (k3_max - dk3) );
end
%2D k-vector grid
[gridObj.k2rep2,gridObj.k1rep2] = meshgrid(gridObj.k2,gridObj.k1);
%Grid indices
gridObj.k1ind0 = floor( n1 ./ 2 ) + 1;
gridObj.k2ind0 = floor( n2 ./ 2 ) + 1;
gridObj.k3ind0 = floor( n3 ./ 2 ) + 1;
end
