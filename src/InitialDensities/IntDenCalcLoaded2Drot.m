%IntDenCalcLoaded2Drot: Initial density is loaded from a saved rho. Fit it to
% you box using linear interpolation

function [ rho ] = IntDenCalcLoaded2Drot( filename, systemObj)
% paths
path = './src/InitialDensities/SavedRhos/';
loadname = [path filename];
% load it. rho should be the name of the saved variable
load( loadname );
% Now fit it to the box 
[NxiTemp, NyiTemp, NmiTemp] = size( rho );
% Put it in a box of size one and wrap rho b/c of PBC
% original size
xiTemp = 0: 1 / NxiTemp: 1 ;
yiTemp = 0: 1 / NyiTemp: 1 ;
miTemp = 0: 1 / NmiTemp: 1 ;
% new size
xfTemp = 0: 1 / systemObj.Nx : 1 ;
yfTemp = 0: 1 / systemObj.Ny : 1 ;
mfTemp = 0: 1 / systemObj.Nm : 1 ;
% build a the 3D variables
[y3fTemp, x3fTemp, m3fTemp] = meshgrid( yfTemp, xfTemp, mfTemp );
% wrap
indsX =  [ 1:NxiTemp 1];
indsY =  [ 1:NyiTemp 1];
indsM =  [ 1:NmiTemp 1];
rhoWrap = rho( indsX, indsY, indsM );
% use n-d linear interp. Don't need 3d variables for inputs, but you do for
% final grid
rho = interpn( xiTemp, yiTemp, miTemp, rhoWrap, x3fTemp, y3fTemp, m3fTemp );
rho = rho( 1:systemObj.Nx, 1:systemObj.Ny, 1:systemObj.Nm );
end

