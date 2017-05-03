%% rotation test
% parameters
Lbox = 1;
n1 = 64;
n2 = 64;
n3 = 8;
sigX = 0.25 * 0.125;
sigY = 0.25 * 0.25;
% grid stuff
% x
dx1 = Lbox / n1;
xInds = 0:n1-1;
xWrapInd = [0:n1/2 -n1/2+1:-1]; 
xWrap = dx1 * xWrapInd;
xCenterInd = -n1/2+1:n1/2;
xCenter = dx1 .* xCenterInd;
% y
dx2 = Lbox / n2;
yInds = 0:n2-1;
yWrapInd = [0:n2/2 -n2/2+1:-1]; 
yWrap = dx2 * yWrapInd;
yCenterInd = -n2/2+1:n2/2;
yCenter = dx2 .* yCenterInd;
% meshgrids
[y2d, x2d] = meshgrid(yWrap,xWrap);
[yCenter2d, xCenter2d] = meshgrid(yCenter,xCenter);
% phi
dPhi = 2*pi / n3;
phi  = 0:dPhi:2*pi-dPhi;
% trig functions
cosPhi = cos(phi);
sinPhi = sin(phi);
% all combinations
indCombWrap = combvec( xWrap, yWrap );
indCombCenter = combvec( xCenter, yCenter );
% allocate
rotIndsWrap = zeros( n1*n2, n3 );
rotIndsCenter = zeros( n1*n2, n3 );
fRotWrap =  zeros( n1, n2, n3 );
fRotCenter =  zeros( n1, n2, n3 );
% funciton to rotate
fWrap = exp( - ( x2d ) .^2 ./ ( 2 * sigX ^ 2 )  - ( y2d ) .^2 ./ ( 2 * sigY ^ 2 )  );
fCenter = exp( - ( xCenter2d ) .^2 ./ ( 2 * sigX ^ 2 )  - ( yCenter2d ) .^2 ./ ( 2 * sigY ^ 2 )  );
% periodic
fWrapPbc = fWrap( [1:n1 1], [1:n2 1]); % include 0 point
fCenterPbc = fCenter( [n1  1:n1], [n2  1:n2]); % include -n/2 point
% loop over phi
for ii = 1:n3
  %rotate
  rotMat = [ cosPhi(ii) sinPhi(ii); -sinPhi(ii) cosPhi(ii) ];
  rotWrap = rotMat * indCombWrap;
  rotCenter = rotMat * indCombCenter;
  % Wrap
  tempXWrap = rotWrap(1,:) ./ dx1; 
  tempYWrap = rotWrap(2,:) ./ dx2; 
  % Fix boundaries
  tempXWrap( tempXWrap < 0 ) = (n1 + tempXWrap( tempXWrap < 0 ) );
  tempYWrap( tempYWrap < 0 ) = (n2 + tempYWrap( tempYWrap < 0 ) );
  % round and store inds
  tempXWrapInds = mod( round( tempXWrap ), n1 ) + 1;
  tempYWrapInds = mod( round( tempYWrap ), n2 ) + 1;
  tempInds =  sub2ind( [n1 n2], tempXWrapInds, tempYWrapInds );
  rotIndsWrap(:,ii) = tempInds;
  % make it a grid
  tempXWrap = reshape( tempXWrap, [n1 n2] );
  tempYWrap = reshape( tempYWrap, [n1 n2] );
  % Center
  tempXCenter = rotCenter(1,:) ./ dx1; 
  tempYCenter = rotCenter(2,:) ./ dx2; 
  % fix wrapping
  tempXCenter( tempXCenter < -n1/2 + 1 ) = (n1 + tempXCenter( tempXCenter < -n1/2 + 1 ) );
  tempXCenter( tempXCenter > n1/2 ) = ( tempXCenter( tempXCenter > n1/2 )  - n1 );
  tempYCenter( tempYCenter < -n2/2 + 1 ) = (n2 + tempYCenter( tempYCenter < -n2/2 + 1 ) );
  tempYCenter( tempYCenter > n2/2 ) = ( tempYCenter( tempYCenter > n2/2 )  - n2 );
  % round and store inds
  tempXCenterInds = mod( round( tempXCenter ) + n1/2 -1, n1 ) + 1;
  tempYCenterInds = mod( round( tempYCenter ) + n2/2 -1, n2 ) + 1;
  tempInds =  sub2ind( [n1 n2], tempXCenterInds, tempYCenterInds );
  rotIndsCenter(:,ii) = tempInds;
  % make it a grid
  tempXCenter = reshape( tempXCenter, [n1 n2] );
  tempYCenter = reshape( tempYCenter, [n1 n2] );
  % Wrap
  fTemp = interpn( [xInds n1],  [yInds n2], ...
    fWrapPbc, tempXWrap, tempYWrap );  
  fRotWrap(:,:,ii) = fTemp(1:end, 1:end);
  % center
  fTemp = interpn( [-n1/2 xCenterInd], [-n2/2 yCenterInd], ...
    fCenterPbc, tempXCenter, tempYCenter );  
  fRotCenter(:,:,ii) = fTemp(1:end, 1:end);
end
% plot centered interp
figure()
for ii = 1:n3
  subplot(2,4,ii);
  imagesc( fRotCenter(:,:,ii) );
  title( ['interp center phi = ' num2str( phi(ii) ) ] );
end
% plot wrap interp
figure()
for ii = 1:n3
  subplot(2,4,ii)
  imagesc( fRotWrap(:,:,ii) );
  title( ['interp wrap phi = ' num2str( phi(ii) ) ] );
end
% plot centered inds
figure()
for ii = 1:n3
  subplot(2,4,ii);
  imagesc( fCenter( reshape( rotIndsCenter(:,ii), [n1 n2] ) ) );
  title( ['inds center phi = ' num2str( phi(ii) ) ] );
end
% plot wrap inds
figure()
for ii = 1:n3
  subplot(2,4,ii)
  imagesc( fWrap( reshape( rotIndsWrap(:,ii), [n1 n2] ) ) );
  title( ['inds wrap phi = ' num2str( phi(ii) ) ] );
end