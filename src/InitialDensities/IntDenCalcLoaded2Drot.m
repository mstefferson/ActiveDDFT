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

xiTemp = 0: 1 / NxiTemp: (1 - 1 / NxiTemp );
yiTemp = 0: 1 / NyiTemp: (1 - 1 / NyiTemp );
miTemp = 0: 1 / NmiTemp: (1 - 1 / NmiTemp );

xfTemp = 0: 1 / systemObj.Nx : (1 - 1 / systemObj.Nx );
yfTemp = 0: 1 / systemObj.Ny : (1 - 1 / systemObj.Ny );
mfTemp = 0: 1 / systemObj.Nm : (1 - 1 / systemObj.Nm );

% build a the 3D variables
[y3fTemp, x3fTemp, m3fTemp] = meshgrid( yfTemp, xfTemp, mfTemp );

% use n-d linear interp. Don't need 3d variables for inputs, but you do for
% final grid
rho = interpn( xiTemp, yiTemp, miTemp, rho, x3fTemp, y3fTemp, m3fTemp );

end

