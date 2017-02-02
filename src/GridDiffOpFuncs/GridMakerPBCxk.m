function [gridObj] = GridMakerPBCxk(n1,n2,n3,l1,l2,l3)

% Make pos grid spacings
dx   = l1/n1;
dy   = l2/n2;
dphi = l3/n3;

% Make vectors and grids
gridObj.x1 = ( 0: dx : l1 - dx );
gridObj.x2 = ( 0: dy : l2 - dy );
gridObj.x3 = ( 0: dphi: (l3 - dphi) );

% Make k-space spacings
dkx = 2*pi/l1;
dky = 2*pi/l2;
dkm = 2*pi/l3; %  = 1

% Make k vectors and grids
kx_max = pi/dx; %Maximum spatial k-vector allowed by grid
ky_max = pi/dy; %Maximum spatial k-vector allowed by grid
km_max = pi / dphi; %Maximum angular k-vector

% For generalization purposes, handle N = 0 cases differently.
if n1 == 1
  gridObj.k1 = 0;
else
  gridObj.k1 = ( -kx_max: dkx: (kx_max - dkx) );
end
if n2 == 1
  gridObj.k2 = 0;
else
  gridObj.k2 =  ( -ky_max: dky: (ky_max - dky) );
end
if n3 == 1
  gridObj.km = 0;
else
  gridObj.km = ( -km_max: dkm: (km_max - dkm) );
end

%2D k-vector grid
[gridObj.k2rep2,gridObj.k1rep2] = meshgrid(gridObj.k2,gridObj.k1);

end
