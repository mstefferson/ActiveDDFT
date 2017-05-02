% %%
% A = rand(64,64,64);
% 
% tic 
% for ii = 1:64
% Aft = fftshift( fftn(A) );
% end
% toc

%% rotation test
Lbox = 1;
n1 = 64;
n2 = 64;
n3 = 8;
sigX = 0.25 * 0.125;
sigY = 0.25 * 0.25;
dx1 = Lbox / n1;
x = dx1  .* [0:n1/2 -n1/2+1:-1];
xIndCenter = [-n1/2+1:n1/2 ];
xCenter = dx1 .* xIndCenter;
% x = -Lbox/2:dx1:Lbox/2 -dx1;
dx2 = Lbox / n2;
y = dx2  .* [0:n2/2 -n2/2+1:-1];
yIndCenter = [-n1/2+1:n1/2 ];
yCenter = dx1 .* yIndCenter;

[y2, x2] = meshgrid(y,x);
[yCenter2, xCenter2] = meshgrid(yCenter,xCenter);
x1Inds = 1:n1;
x2Inds = 1:n2;
x3Inds = 1:n3;

dPhi = 2*pi / n3;
phi  = 0:dPhi:2*pi-dPhi;

cosPhi = cos(phi);
sinPhi = sin(phi);
%% Rotate one point
% xTemp = 0;
% yTemp = 6;
% singRot = zeros( 1, n3 );
% for ii = 1:n3
%   rotMat = [ cosPhi(ii) sinPhi(ii); -sinPhi(ii) cosPhi(ii) ];
%   rotVecs = round( rotMat * [xTemp; yTemp] );
%   % keep zero inds
%   tempX = mod( rotVecs(1,:)', n1 ); 
%   tempY = mod( rotVecs(2,:)', n2 );
%   % shift back
%   smartInd = sub2ind( [n1 n2], tempX+1, tempY+1 );
%   singRot(ii) = smartInd;
%   matTemp = zeros(n1,n2);
%   matTemp(smartInd) = 1;
%   subplot(2,4,ii)
%   imagesc(matTemp)
%   axis square
%   title( ['phi = ' num2str( phi(ii) ) ] );
% end

%% Roate all points
% deltaVec = combvec( 0:n1-1, 0:n2-1);
% periodic method

xDeltaWrap = [0:n1/2 -n1/2+1:-1]; 
xDelta = [0:n1-1]; 
xDeltaPbc = [1:n1 1];
%xDeltaIndsPbc = [-n1/2:-1 0:n1/2 ]; 
%xIndsPbc = [1:n1/2+1 n1/2+1:n1];
% yDeltaPbc = [0:n2/2 -n2/2+1:-1];
yDeltaPbc = [1:n2 1];
yDeltaWrap = [0:n2/2 -n2/2+1:-1]; 
yDelta = [0:n2-1]; 
%IndsPbc = [1:n2/2+1 n2/2+1:n2];
deltaWrapVec = combvec( xDeltaWrap, yDeltaWrap );
deltaCenter = combvec( xIndCenter, yIndCenter );
deltaVec = combvec( xDeltaPbc, yDeltaPbc );
% rep
%[yDelta2, xDelta2] = meshgrid(yDeltaIndsPbc,xDeltaIndsPbc);
% deltaXRot = zeros( n1, n2, n3 );
% deltaYRot = zeros( n1, n2, n3 );
newInds = zeros( n1*n2, n3 );
fNew =  zeros( n1, n2, n3 );
fCenterNew =  zeros( n1, n2, n3 );
f = exp( - ( x2 ) .^2 ./ ( 2 * sigX ^ 2 )  - ( y2 ) .^2 ./ ( 2 * sigY ^ 2 )  );
fCenter = exp( - ( xCenter2 ) .^2 ./ ( 2 * sigX ^ 2 )  - ( yCenter2 ) .^2 ./ ( 2 * sigY ^ 2 )  );

fPbc = f(xDeltaPbc, yDeltaPbc);
fCenterPbc = fCenter( [n1  1:n1], [n2  1:n2]);

for ii = 1:n3
  rotMat = [ cosPhi(ii) sinPhi(ii); -sinPhi(ii) cosPhi(ii) ];
%   rotVecs = round( rotMat * deltaVec2 );
  %rotNoRound = rotMat * deltaVecPbc;
  rotIndsNoRound = rotMat * deltaWrapVec;
  rotIndsCenter = rotMat * deltaCenter;
%   xRot = mod( reshape(rotVecs(1,:), [n1 n2] ), n1);
%   yRot = mod( reshape(rotVecs(2,:), [n1 n2] ), n2 );
%   tempX = mod( rotVecs(1,:)', n1-1); 
%   tempY = mod( rotVecs(2,:)', n2-1 );
%   tempXnoRound = mod( rotNoRound(1,:)', n1-1 ); 
%   tempYnoRound = mod( rotNoRound(2,:)', n2-1 );
  tempXnoRound = rotIndsNoRound(1,:); 
  tempYnoRound = rotIndsNoRound(2,:); 
  tempXCenter = rotIndsCenter(1,:); 
  tempYCenter = rotIndsCenter(2,:); 
  % check these
  %tempXnoRound( tempXnoRound > n1/2 ) = -(n1 - tempXnoRound( tempXnoRound > n1/2 ) );
  %tempXnoRound( tempXnoRound < -n1/2 ) = (n1 + tempXnoRound( tempXnoRound < -n1/2 ) );

  %tempYnoRound( tempYnoRound > n2/2 ) = tempYnoRound( tempYnoRound > n2/2 ) - n2 + 1;
  %tempYnoRound( tempYnoRound < -n2/2 ) = tempYnoRound( tempYnoRound < -n2/2 ) + n2 -1;
  %
  tempXnoRound( tempXnoRound < 0 ) = (n1 + tempXnoRound( tempXnoRound < 0 ) );
  tempYnoRound( tempYnoRound < 0 ) = (n2 + tempYnoRound( tempYnoRound < 0 ) );
  %
  tempXCenter( tempXCenter < -n1/2 + 1 ) = (n1 + tempXCenter( tempXCenter < -n1/2 + 1 ) );
  tempXCenter( tempXCenter > n1/2 ) = ( tempXCenter( tempXCenter > n1/2 )  - n1 );
  tempYCenter( tempYCenter < -n2/2 + 1 ) = (n2 + tempYCenter( tempYCenter < -n2/2 + 1 ) );
  tempYCenter( tempYCenter > n2/2 ) = ( tempYCenter( tempYCenter > n2/2 )  - n2 );
  %
  xNoRound = reshape( tempXnoRound, [n1 n2] );
  yNoRound = reshape( tempYnoRound, [n1 n2] );
  %
  xCenterTemp = reshape( tempXCenter, [n1 n2] );
  yCenterTemp = reshape( tempYCenter, [n1 n2] );
%   temp = [tempX'; tempY'];
%   smartInd = sub2ind( [n1 n2], tempX+1, tempY+1 );
%   
%   xI = [ 0:n1-1 ];
%  xIwrap = [ 1:n1 1];
%   yI = [ 0:n2-1  ];
%   yIwrap = [ 1:n2 1];
%   fWrap = f(xIwrap,yIwrap);
%   fTemp = interpn( xI, yI, fWrap, xNoRound, yNoRound);

% Do everything in terms of inds
 
  fTemp = interpn( dx1.* [xDelta n1], dx2 .* [yDelta n2], ...
    fPbc, dx1.* xNoRound, dx2 .* yNoRound );  
  fNew(:,:,ii) = fTemp(1:end, 1:end);
  fTemp = interpn( dx1.* [-n1/2 xIndCenter], dx2 .* [-n2/2 yIndCenter], ...
    fCenterPbc, dx1.* xCenterTemp, dx2 .* yCenterTemp );  
  fCenterNew(:,:,ii) = fTemp(1:end, 1:end);
end

figure()
for ii = 1:n3
  subplot(2,4,ii);
  imagesc( fCenterNew(:,:,ii) );
  title( ['phi = ' num2str( phi(ii) ) ] );
end

figure()
for ii = 1:n3
  subplot(2,4,ii)
  imagesc( fNew(:,:,ii) );
  title( ['phi = ' num2str( phi(ii) ) ] );
end



