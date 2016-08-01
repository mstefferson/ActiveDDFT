function [gridObj] = GridMakerPBCxk(Nx,Ny,Nm,Lx,Ly,Lphi)

% Make pos grid spacings
dx   = Lx/Nx;              
dy   = Ly/Ny;              
dphi = Lphi/Nm; 

% Make vectors and grids
gridObj.x = ( 0: dx : Lx - dx );
gridObj.y = ( 0: dy : Ly - dy );
gridObj.phi = ( 0: dphi: (Lphi - dphi) );

% Make k-space spacings
dkx = 2*pi/Lx;
dky = 2*pi/Ly;
dkm = 2*pi/Lphi; %  = 1

% Make k vectors and grids
kx_max = pi/dx; %Maximum spatial k-vector allowed by grid
ky_max = pi/dy; %Maximum spatial k-vector allowed by grid
km_max = pi / dphi; %Maximum angular k-vector

% For generalization purposes, handle N = 0 cases differently.
if Nx == 0
  grid.kx = 0;
else
  gridObj.kx = ( -kx_max: dkx: (kx_max - dkx) );
end
if Ny == 0
  gridObj.ky = 0;
else
  gridObj.ky =  ( -ky_max: dky: (ky_max - dky) );
end
if Nm == 0 
  gridObj.km = 0;
else
gridObj.km = ( -km_max: dkm: (km_max - dkm) );
end

%2D k-vector grid
[gridObj.ky2D,gridObj.kx2D] = meshgrid(gridObj.ky,gridObj.kx);

end
