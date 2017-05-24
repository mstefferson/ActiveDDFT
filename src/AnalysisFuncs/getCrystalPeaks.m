% [peaks] = getCrystalPeaks( rho, lBox, plotFlag )
%
% Get the k-space peaks of a concentration
%
function [peaks] = getCrystalPeaks( rho, lBox, plotFlag )
% make sure size is ok
[n1, n2]  = size( rho );
if ( n1 > 1 ) && ( n1 > 1 ) && ( n1 ~= n2 )
  error('Only working for square boxes');
end
n = max( n1, n2 );
dk = 2 * pi / lBox;
if n > 1
  numPeaksK = min( 7, n ); % max peaks to look for
  % positive k vectors
  k1 = 0:dk:pi*n1/lBox-dk;
  k2 = 0:dk:pi*n2/lBox-dk;
  % FT
  finalFT = fftshift( fftn( rho ) );
  centerInd1 = floor( n1 / 2 ) + 1;
  centerInd2 = floor( n2 / 2 ) + 1;
  % number of mode cutoff
  numModes1 = min( 12, n1-1);
  numModes2 = min( 12, n2-1);
  inds1 = centerInd1 - numModes1 : centerInd1 + numModes1;
  inds2 = centerInd2 - numModes2 : centerInd2 + numModes2;
  numInds1 = length(inds1);
  numInds2 = length(inds2);
  desiredFT = finalFT(inds1,inds2);
  desiredFTSq = abs( desiredFT ) .^ 2;
  % grab peaks
  peakFrac2Keep = 0.01;
  [topPeaks, topPeakInds] = sort( desiredFTSq(:), 'descend');
  % get average ignoring first peak
  topPeakInds = topPeakInds(1:numPeaksK);
  topPeakInds = topPeakInds( topPeaks(1:numPeaksK) > peakFrac2Keep * topPeaks(1) );
  numPeaksK = length( topPeakInds );
  [iP, jP] = ind2sub( [ numInds1 numInds2 ], topPeakInds');
  peakValK = topPeaks(1:numPeaksK);
  %  get peak k values and calc distance
  kVecAll = [ ( k1(iP(1)) - k1(iP) )' ( k2(jP(1)) - k2(jP) )' ];
  dist2PeaksAll = sqrt( kVecAll(:,1).^2 + kVecAll(:,2).^2 );
  indsKAll = [ (iP(1) - iP)' (jP(1) - jP)' ];
  % just positive
  conds1 =  k1( iP(1) ) - k1(iP) >= 0;
  conds2 = k2( jP(1) ) - k2(jP) >= 0;
  conds = conds1 .* conds2;
  iP = iP( conds == 1 );
  jP = jP( conds == 1 );
  kVecPos = [ ( k1(iP(1)) - k1(iP) )' ( k2(jP(1)) - k2(jP) )' ];
  dist2PeaksPos = sqrt( kVecPos(:,1).^2 + kVecPos(:,2).^2 );
  indsKPos = [ (iP(1) - iP)' (jP(1) - jP)' ];
  ka = mean( dist2PeaksPos( dist2PeaksPos > 0 ) );
  % real space lattice spacing
  if numPeaksK > 6
    a =  4 * pi / ( sqrt(3) * mean( unique( dist2PeaksPos( dist2PeaksPos > 0) ) ) ) ;
  elseif numPeaksK == 3
    a = 2 * pi / mean( unique( dist2PeaksPos( dist2PeaksPos > 0 ) ) );
  else
    a = NaN;
  end
  % plotting
  if plotFlag
    % get real position peaks
    if ~isnan( a ) && a < lBox / 2
      % grab some section in the middle
      deltaX = ceil( a / lBox * n1 );
      inds1 = centerInd1 - deltaX :centerInd1 + deltaX ;
      inds2 = centerInd2 - deltaX :centerInd2 + deltaX ;
      desiredReal = rho(inds1, inds2);
      newSize = size(desiredReal);
      [~,newMax] = max( desiredReal(:) );
      [iMax, jMax] = ind2sub( newSize , newMax ) ;
      % Center it on a peak
      newCenter1 = inds1(iMax);
      newCenter2 = inds2(jMax);
      inds1 = newCenter1 - deltaX : newCenter1 + deltaX;
      inds2 = newCenter2 - deltaX : newCenter2 + deltaX;
      desiredReal = rho(inds1, inds2);
    else
      desiredReal = rho;
    end
    % handle grid stuff
    dx = lBox / n;
    sizeReal = size( desiredReal );
    if min(n1,n2) > 1
      xVec = dx * (-ceil( ( max(sizeReal) - 1 ) / 2 ): floor( ( max(sizeReal) - 1 ) / 2 ) );
    else
      xVec = dx * (0:n-1);
    end
    sizeK = size( desiredFTSq );
    kVec = dk * (-ceil( ( max(sizeK) - 1 ) / 2 ): floor( ( max(sizeK) - 1 ) / 2 ) );
    % plot time
    figure()
    if min(n1,n2) > 1
      subplot(1,2,1);
      imagesc( xVec, xVec, desiredReal' );
      xlabel('$$ x $$'); ylabel('$$ y $$');
      title( 'real');
      colorbar
      subplot(1,2,2);
      imagesc( kVec, kVec, desiredFTSq' );
      xlabel('$$ k_x $$'); ylabel('$$ k_y $$');
      title( 'k-space');
      colorbar
    else
      subplot(1,2,1);
      plot( xVec, desiredReal );
      xlabel('$$ x $$'); ylabel('$$ f $$');
      title( 'real');
      subplot(1,2,2);
      plot( kVec, desiredFTSq );
      xlabel('$$ k_x $$'); ylabel('$$ | f(k) | ^ 2  $$');
      title( 'k-space');
    end
  end
else
  % nothing to find
  numPeaksK = 1;
  peakValK = rho;
  kVecAll = [0 0];
  dist2PeaksAll = 0;
  indsKAll = [0 0];
  kVecPos = [0 0];
  dist2PeaksPos = 0;
  indsKPos = [0 0];
  a = 0;
  ka = 0;
end
% Store everything
peaks.numPeaksK = numPeaksK;
peaks.peakValK = peakValK;
peaks.kVecAll = kVecAll;
peaks.dist2PeaksAll = dist2PeaksAll;
peaks.indsKAll = indsKAll;
peaks.kVecPos = kVecPos;
peaks.dist2PeaksPos = dist2PeaksPos';
peaks.indsKPos = indsKPos;
peaks.dk = dk;
peaks.a = a;
peaks.ka = ka;
