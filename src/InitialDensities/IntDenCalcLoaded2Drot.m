%IntDenCalcLoaded2Drot: Initial density is loaded from a saved rho. Fit it to
% you box using linear interpolation

function [ rho ] = IntDenCalcLoaded2Drot( filename, systemObj)
% paths
path2file = './src/InitialDensities/SavedRhos/';
% path2file = './';
loadname = [path2file filename];
% load it. rho should be the name of the saved variable
load( loadname );
% Now fit it to the box 
[n1iTemp, n2iTemp, n3iTemp] = size( rho );
% Put it in a box of size one and wrap rho b/c of PBC
% original size
xiTemp = 0: 1 / n1iTemp: 1 ;
yiTemp = 0: 1 / n2iTemp: 1 ;
miTemp = 0: 1 / n3iTemp: 1 ;
% new size
xfTemp = 0: 1 / systemObj.n1 : 1 ;
yfTemp = 0: 1 / systemObj.n2 : 1 ;
mfTemp = 0: 1 / systemObj.n3 : 1 ;
% build a the 3D variables
[y3fTemp, x3fTemp, m3fTemp] = meshgrid( yfTemp, xfTemp, mfTemp );
% wrap
indsX =  [ 1:n1iTemp 1];
indsY =  [ 1:n2iTemp 1];
indsM =  [ 1:n3iTemp 1];
rhoWrap = rho( indsX, indsY, indsM );
% use n-d linear interp. Don't need 3d variables for inputs, but you do for
% final grid
rho = interpn( xiTemp, yiTemp, miTemp, rhoWrap, x3fTemp, y3fTemp, m3fTemp );
rho = rho( 1:systemObj.n1, 1:systemObj.n2, 1:systemObj.n3 );
end

