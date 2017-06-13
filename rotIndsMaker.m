
function [angleIndMap, rotInds] = rotIndsMaker( n1, n2, n3, l1, l2, lRod )

% x1
dx1 = l1 / n1;
x1Inds = 0:n1-1;
x1WrapInd = [0:n1/2 -n1/2+1:-1];
x1Wrap = dx1 * xWrapInd;
% x2
dx2 = l2 / n2;
x2Inds = 0:n2-1;
x2WrapInd = [0:n2/2 -n2/2+1:-1];
x2Wrap = dx2 * x2WrapInd;
indCutoff1 = ceil( lRod / l1 ) .* n1  + 1;
indCutoff2 = ceil( lRod / l2 ) .* n2  + 1;
if indCutoff1 > n1/2 || indCutoff2 > n2/2
  error( 'lrod / l1  is too large to do Mayer function rotation' )
  fprintf( 'lrod / l1  is too large to do Mayer function rotation \n' )
end
% get sub indices
inds2rotWrap1 = [ 1:(indCutoff1+1) (n1-indCutoff1+1):n1];
inds2rotWrap2 = [ 1:(indCutoff2+1) (n2-indCutoff2+1):n2];
combRodInds = combvec( inds2rotWrap1, inds2rotWrap2 ); 
inds2rotList = sub2ind( [n1 n2], combRodInds(1,:)', combRodInds(2,:)' ); 
% phi
dPhi = 2*pi / n3;
phi  = 0:dPhi:2*pi-dPhi;
% trig functions
cosPhi = cos( phi(1:n3/2) );
sinPhi = sin( phi(1:n3/2) );
% all combinations
indCombWrap = combvec( xWrap, yWrap );
indCombCenter = combvec( xCenter, yCenter );
% reused indices for mapping
numRot = n3/2;
angleIndMap = [1:numRot 1:numRot];
unrotInds = (1:n1*n2)';
rot3rep = reshape( repmat( 1:n3, [n1*n2,1] ), [n1*n2*n3, 1] );
% allocate
rotInds = zeros(n1*n2*n3,n3);
% loop over phi
for ii = 1:numRot
  %rotate
  rotMat = [ cosPhi(ii) sinPhi(ii); -sinPhi(ii) cosPhi(ii) ];
  rotWrap = rotMat * indCombWrap;
  % Wrap
  tempXWrap = rotWrap(1,:) ./ dx1;
  tempYWrap = rotWrap(2,:) ./ dx2;
  % Fix boundaries
  tempXWrap( tempXWrap < 0 ) = (n1 + tempXWrap( tempXWrap < 0 ) );
  tempYWrap( tempYWrap < 0 ) = (n2 + tempYWrap( tempYWrap < 0 ) );
  % round and store inds
  tempXWrapInds = round( tempXWrap ) + 1;
  tempYWrapInds = round( tempYWrap ) + 1;
  tempInds =  sub2ind( [n1 n2], tempXWrapInds, tempYWrapInds );
  % subset
  % smart indices
  newInds = unrotInds;
  newInds(inds2rotList) = tempInds(inds2rotList);
  % in x,y form
  [rot1, rot2] =  ind2sub( [n1 n2], newRotInds );
  rot1rep = repmat( rot1, [n3, 1] );
  rot2rep = repmat( rot2, [n3, 1] );
  tempInds =  sub2ind( [n1 n2 n3], rot1rep, rot2rep, rot3rep );
  rotInds(:,ii) = tempInds;
end
