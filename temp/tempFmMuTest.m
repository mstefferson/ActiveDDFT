% params
n = 64;
n1 = n;
n2 = n;
n3 = n;
l1 = 3;
l2 = 3;
lRod = 1;
% rho
rho = 3 .* ones( n1, n2, n3 );
rhoFt = fftshift( fftn( rho ) );
% system obj
systemObj.n1 = n1;
systemObj.n2 = n2;
systemObj.n3 = n3;
systemObj.l1 = l1;
systemObj.l2 = l2;
systemObj.l3 = 2*pi;
% inds
inds = 1:n3;
indsMinus = [1 n3:-1:2 ];
centerInd1 = n1/2 + 1;
centerInd2 = n2/2 + 1;
centerInd3 = n3/2 + 1;
% scale factor
scaleFact = (  systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
  (systemObj.n1 * systemObj.n2 * systemObj.n3 ^ 2);
% mayer
[fm, fmFt] = mayerFncHrLabFrame( n1, n2, n3, l1, l2, lRod );
% muEx
[muExFt1] = muExCalcMayerLF(rhoFt,fmFt,systemObj,scaleFact, inds, indsMinus);
muEx = ifftn( ifftshift( muExFt1 ) );

% allocate
% muExFt = zeros( systemObj.n1, systemObj.n2, systemObj.n3 );
% % sum over angular modes
% for ll = 1:systemObj.n3
%   muExFt = muExFt + repmat( rhoFt( :, :, inds(ll) ), [1 1 systemObj.n3] ) .* ...
%     reshape( fmFt( :, :, :, indsMinus(ll) ), [systemObj.n1 systemObj.n2 systemObj.n3] );
% end
% muExFt = -scaleFact .* muExFt;
% look at center
ll = centerInd3;

fmFtCenter = reshape( fmFt( centerInd1, centerInd2, :, centerInd3 ), [1 systemObj.n3] );
fmCenter = ifftn( ifftshift( fmFtCenter ) );
fmFtReshape = reshape( fmFt( :, :, :, indsMinus(ll) ), [systemObj.n1 systemObj.n2 systemObj.n3] );
fmReshape = ifftn( ifftshift( fmFtReshape ) );
