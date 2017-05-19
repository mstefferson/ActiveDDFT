%%
function [peaks] = getCrystalPeaks( rho, 1 );
% make sure size is ok
[n1, n2]  = size( rho, 1 );
if n1 ~= n2
  error('Only written for square boxes');
else
  n = n1;
end
% FT
finalFT = fftshift( fftn( rho ) );
centerInd = n / 2 + 1;
numModes = 12;
inds = centerInd - numModes : centerInd + numModes;
numInds = length(inds);
desiredFT = finalFT(inds,inds);
desiredFTSq = abs( desiredFT ) .^ 2;
%% grab peaks
numPeaks = 7;
meanPeaks = mean( desiredFTSq(:) );
stdPeaks = std( desiredFTSq(:) );
[topPeaks, topPeakInds] = sort( desiredFTSq(:), 'descend');
topPeakInds = topPeakInds(1:numPeaks);
topPeakInds = topPeakInds( topPeaks(1:numPeaks) > meanPeaks + stdPeaks );
numPeaks = length( topPeakInds );
[iP, jP] = ind2sub( [ numInds numInds ], topPeakInds');
% distance calc
dist2Peaks =  sqrt( ( iP(1) - iP ) .^ 2  + ( jP(1) - jP ) .^ 2 );
% get real position peaks
% grab some section in the middle
numSites = 50;
a =  2 * l / ( sqrt(3) * mean( unique( dist2Peaks(2:end) ) ) ) ;
deltaX = ceil( a / l * n );
inds = centerInd - deltaX :centerInd + deltaX ;
desiredReal = denRecObj.rhoFinal(inds, inds);
newSize = size(desiredReal);
[~,newMax] = max( desiredReal(:) );
[iMax, jMax] = ind2sub( newSize , newMax ) ;
% Center it on a peak
newCenterX = inds(iMax);
newCenterY = inds(jMax);
indsX = newCenterX - deltaX : newCenterX + deltaX;
indsY = newCenterY - deltaX : newCenterY + deltaX;
desiredReal = denRecObj.rhoFinal(indsX, indsY);
newCenter = numSites + 1;
% Store everything
peaks.numPeaksK = numPeaks;
peaks.peakValK = topPeaks(1:numPeaks);
peaks.distK = 2 * pi / l * dist2Peaks;
peaks.indsK = [ (iP-iP(1))' (jP-jP(1))'];
peaks.a = a;
peaks.ka = mean( unique( dist2Peaks(2:end) ) ) * 2 * pi / l;
% plotting
% handle grid stuff
dx = l / n;
dk = 2*pi / l;
sizeReal = size( desiredReal, 1  );
dReal = round( ( sizeReal - 1 ) / 2 );
xVec = (-dReal:dReal ) * dx ;
sizeK = size( desiredFTSq, 1  );
dK = round( ( sizeK - 1 ) / 2 );
kVec = ( -dK:dK ) * dk;
% plot time
figure()
subplot(1,2,1);
imagesc( xVec, xVec, desiredReal );
title( 'real');
subplot(1,2,2);
imagesc( kVec, kVec, desiredFTSq );
title( 'k-space');
