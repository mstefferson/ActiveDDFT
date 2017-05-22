%%
function [peaks] = getCrystalPeaks( rho, lBox, plotFlag )
% make sure size is ok
[n1, n2]  = size( rho );
if ( n1 > 1 ) && ( n1 > 1 ) && ( n1 ~= n2 )
  error('Only working for square boxes');
end
n = max( n1, n2 );
if n > 1
  maxNumPeaks = min( 7, n ); % max peaks to look for
  dk = 2 * pi / lBox;
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
  %% grab peaks
  peakFrac2Keep = 0.01;
  [topPeaks, topPeakInds] = sort( desiredFTSq(:), 'descend');
  % get average ignoring first peak
  %   meanPeaks = mean( topPeaks(2:end) );
  %   stdPeaks = std( topPeaks(2:end) );
  topPeakInds = topPeakInds(1:maxNumPeaks);
  topPeakInds = topPeakInds( topPeaks(1:maxNumPeaks) > peakFrac2Keep * topPeaks(1) );
  maxNumPeaks = length( topPeakInds );
  [iP, jP] = ind2sub( [ numInds1 numInds2 ], topPeakInds');
  % distance calc
  dist2Peaks =  dk .* sqrt( ( iP(1) - iP ) .^ 2  + ( jP(1) - jP ) .^ 2 );
  % real space lattice spacing
  if maxNumPeaks > 6
    a =  4 * pi / ( sqrt(3) * mean( unique( dist2Peaks( dist2Peaks > 0) ) ) ) ;
  elseif maxNumPeaks == 3
    a = 2 * pi / mean( unique( dist2Peaks( dist2Peaks > 0 ) ) );
  else
    a = NaN;
  end
  % Store everything
  peaks.numPeaksK = maxNumPeaks;
  peaks.peakValK = topPeaks(1:maxNumPeaks);
  peaks.distK = dist2Peaks;
  peaks.distKunique = unique( dist2Peaks );
  peaks.indsK = [ (iP-iP(1))' (jP-jP(1))' ];
  peaks.dk = dk;
  peaks.a = a;
  peaks.ka = dk * mean( unique( dist2Peaks > 0 ) );
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
  % Store everything
  peaks.numPeaksK = 1;
  peaks.peakValK = rho;
  peaks.distK = 0;
  peaks.distKunique = 0;
  peaks.indsK = 0;
  peaks.dk = 0;
  peaks.a = 0;
  peaks.ka = 0;
end