function [gridObj] = GridMakerPBCxk(Nx,Ny,Nm,Lx,Ly)

% Make pos grid spacings
dx   = Lx/Nx;              
dy   = Ly/Ny;              
dphi = 2*pi/Nm;

% Make vectors and grids
gridObj.x                = ( 0: dx : Lx - dx );
% x                = ( -Lx/2 : dx: Lx/2 - dx); 
gridObj.y                = ( 0: dy : Ly - dy );
% y                = ( -Ly/2 : dy: Ly/2 - dy);
gridObj.phi              = ( 0: dphi: (2*pi - dphi) );

% Make k-space spacings
dkx          = 2*pi/Lx;
dky          = 2*pi/Ly;
dkm          = 1;

% Make k vectors and grids
kx_max           = pi/dx; %Maximum spatial k-vector allowed by grid
ky_max           = pi/dy; %Maximum spatial k-vector allowed by grid
km_max           = pi / dphi; %Maximum angular k-vector
gridObj.kx               = ( -kx_max: dkx: (kx_max - dkx) );
gridObj.ky               = ( -ky_max: dky: (ky_max - dky) );
gridObj.km               = ( -km_max: dkm: (km_max - dkm) );
%2D k-vector grid
[gridObj.ky2D,gridObj.kx2D] = meshgrid(gridObj.ky,gridObj.kx);

end
