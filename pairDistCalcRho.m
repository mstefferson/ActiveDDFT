% [pDist] = pairDistCalcRho(rho,l1,l2,lRod)
% 
% Calculates the density dependent pair correlation function
% g(r,r') = rho^(2) ./ rho^(1)rho^(1)  
%
function [pDistAve] = pairDistCalcRho(rho,l1,l2,lRod)
[n1,n2,n3] = size( rho );
% add paths just in case
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
% build phi
dx1 = l1 ./ n1;
dx2 = l2 ./ n2;
dphi = 2 * pi / n3;
x1 = dx1 .* (-n1/2 + 1:n1/2 ); 
x2 = dx2 .* (-n2/2 + 1:n2/2 ); 
phi = 0 : dphi : 2*pi - dphi;
% get MayerFunction
[mayer] = mayerFncHr(n1, n2, n3, l1, l2, lRod) ;
% calc concentration
c = trapz_periodic( phi, rho, 3 );
% calculate max indices to check
maxDelta1 =  ceil( n1 .* lRod ./ l1 );
indsDelta1 = mod( -maxDelta1:(maxDelta1+1) - 1, n1 ) + 1;
totalInds1 = 2 .* maxDelta1 + 1;
maxDelta2 =  ceil( n2 .* lRod ./ l2 );
indsDelta2 = mod( -maxDelta2:(maxDelta2+1) - 1, n2 ) + 1;
totalInds2 = 2 .* maxDelta2 + 1;
indsDelta3 = 1:n3;
% initialize
pDist = ones( n1, n2,n1, n2);
vec2Int = zeros( totalInds1, totalInds2, n3 );
for ii = 1:n1
  % inds for loop over delta x
  inds1 = mod( (ii-1) + indsDelta1 - 1, n1 ) + 1;
  for jj = 1:n2
    % inds for loop over delta y
    inds2 = mod( (jj-1) + indsDelta2 - 1, n2 ) + 1;
    for mm = 1:n3
      inds3 = mod( (mm-1) + indsDelta3 - 1, n3 ) + 1;
      mat2IntTemp = rho(ii, jj, mm) .* rho( inds1, inds2, inds3)...
        .* mayer( indsDelta1, indsDelta2, indsDelta3 );
      vec2Int(:,:,mm) = trapz_periodic( phi, mat2IntTemp, 3 );
    end
    num = trapz_periodic( phi, vec2Int, 3 );
    den = ( c(ii,jj) .* c(inds1, inds2) );
    newterm = num ./ den;
    reshapeNewterm = reshape( newterm, [1, 1, totalInds1, totalInds2] );
    pDist( ii, jj, inds1, inds2 ) = pDist( ii, jj, inds1, inds2 ) + reshapeNewterm;
  end
end
% Calculate average
pDistAve = ones( n1, n2);
combInd = combvec( 1:n1, 1:n2 );
allInds1 = combInd(1,:).';
allInds2 = combInd(2,:).';
matSize = [n1, n2, n1, n2];
for ii = indsDelta1
  ii2 = mod( allInds1 + (ii-2), n1 ) + 1 ;
  for jj = indsDelta2
    jj2 = mod( allInds2 + (jj-2), n2 ) + 1 ;
    inds = sub2ind( matSize, allInds1', allInds2', ii2',  jj2' );
    tempMean =  mean( pDist( inds ) );
    pDistAve(ii,jj) =  tempMean;
  end
end
% center
center1 = round( n1 / 2 ) + 1 ;
center2 = round( n2 / 2 ) + 1 ;
pDistAve = circshift( circshift( pDistAve, center1-1, 1 ), center2-1, 2 );
% plot
figure()
imagesc( x1, x2, pDistAve )
title('pair distribution g(x,y): density dependent')
xlabel('x'); ylabel('y')
colorbar
keyboard
