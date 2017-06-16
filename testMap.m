clear
clc
% stretch amount
a = 1/2;
%% 1D
% grid stuff
ni = 256;
l = 1;
dk  = 2 * pi / l;
k  = 4 * dk;
dxi = l / ni;
xMax = l / 2;
xi  = -xMax:dxi:xMax-dxi;
xiCenterInd = ni/2+1;
nPad = max( 0, ( ni * ( l / (2*a) - 1 ) ) / 2 );
% new grid
nf = ni + 2*nPad;
xfCenterInd = nf/2+1;
lf = l / a;
dxf = lf / nf;
xf = -lf/2:dxf:lf/2-dxf;
% create function
sigma = l / 10;
f = exp( - xi .^ 2 ./ ( 2 * sigma ^ 2 ) ) .* cos( k * xi );
fPad = zeros( 1, nf );
fPad( 1+nPad:nPad + ni ) = f;
% linear interp
fMap = interp1( xf, fPad, xi );
fMap( isnan( fMap ) ) = 0;
% plot
figure()
plot(xi,f);
hold all
plot(xi,fMap);
legend('original', 'stretched')
%
%% 2D
% grid stuff
k1  = 4 * dk;
k2  = 8 * dk;
[xi2Mesh, xi1Mesh] = meshgrid( xi, xi );
% new grid
[xf2Mesh, xf1Mesh] = meshgrid( xf, xf );
% create function
sigma1 = l / 10;
sigma2 = l / 12;
f = exp( -xi1Mesh .^ 2 ./ ( 2 * sigma1 ^ 2 ) - xi2Mesh .^ 2 ./ ( 2 * sigma2 ^ 2 ) ) .* ...
  cos( k1 * xi1Mesh ) .* cos( k2 * xi2Mesh );
% try rhoFt
path2file = '/Users/mike/Projects/HardRodML/analyzedfiles/analyzeMe/l10/modes16/Hr_disks_softshoulder_1il1_2il1.85_1ie1_2ie2.5_diag1_N2562561_ls1010_bc3.14_vD0_IC0_SM6_t06.01/';
temp = 'params_Hr_disks_softshoulder_1il1_2il1.85_1ie1_2ie2.5_diag1_N2562561_ls1010_bc3.14_vD0_IC0_SM6_t06.01.mat';
load( [path2file temp] );
rhoFt = fftshift( fftn( denRecObj.rhoFinal ) );
f = rhoFt;
%
fPad = zeros( nf, nf );
fPad( 1+nPad:nPad + ni, 1+nPad:nPad + ni ) = f;
% linear interp
fMap = interp2( xf2Mesh, xf1Mesh, fPad, xi2Mesh, xi1Mesh );
fMap( isnan( fMap ) ) = 0;
% plot
figure()
subplot(2,2,1)
imagesc(abs( f ) .^ 2);
colorbar
title('original');
subplot(2,2,2)
imagesc( abs( fMap ) .^ 2 );
colorbar
title('stretched');
subplot(2,2,3)
plot(xi, abs( f(xiCenterInd,:) ) .^ 2 , xi, abs( f(:,xiCenterInd) ) .^ 2)
legend('columns','rows','location','best')
title('original')
subplot(2,2,4)
plot(xi, abs( fMap(xiCenterInd,:) ) .^2 , xi, abs( fMap(:,xiCenterInd) .^ 2 ) )
title('original')
legend('columns','rows','location','best')
