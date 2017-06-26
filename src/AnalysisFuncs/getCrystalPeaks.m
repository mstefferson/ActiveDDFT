% [peaks] = getCrystalPeaks( rho, lBox, plotFlag )
%
% Get the k-space peaks of a concentration
%
function [cryPeaks] = getCrystalPeaks( rho, lBox, plotFlag )
try
  % make sure size is ok
  [n1, n2]  = size( rho );
  if ( n1 > 1 ) && ( n1 > 1 ) && ( n1 ~= n2 )
    error('Only working for square boxes');
  end
  n = max( n1, n2 );
  dk = 2 * pi / lBox;
  if n > 1
    numPeaksK = min( 7, n ); % max peaks to look for
    % FT
    finalFT = fftshift( fftn( rho ) );
    centerInd1 = floor( n1 / 2 ) + 1;
    centerInd2 = floor( n2 / 2 ) + 1;
    % number of mode cutoff
    % crystal b has smalled spacing a = 1.21. Convert that to mode number
    aMin = 1.21;
    nMode = ceil( 2 * lBox / (sqrt(3)*aMin) );
    numModes1 = min( 2*nMode, n1-1);
    numModes2 = min( 2*nMode, n2-1);
    inds1 = centerInd1 - numModes1 : centerInd1 + numModes1;
    inds2 = centerInd2 - numModes2 : centerInd2 + numModes2;
    % k1 k2 for these inds
    k1Desired =  dk .* (inds1-centerInd1);
    k2Desired =  dk .* (inds2-centerInd2);
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
    kVecAll = [ ( k1Desired(iP) - k1Desired(iP(1)) )' ...
      ( k2Desired(jP) - k2Desired(jP(1)) )' ];
    dist2PeaksAll = sqrt( kVecAll(:,1).^2 + kVecAll(:,2).^2 );
    indsKAll = [ (iP-iP(1))' (jP-jP(1))' ];
    % just positive
    conds1 =  k1Desired(iP) - k1Desired( iP(1) ) >= 0;
    conds2 = k2Desired(jP) - k2Desired( jP(1) )  >= 0;
    conds = conds1 .* conds2;
    iP = iP( conds == 1 );
    jP = jP( conds == 1 );
    kVecPos = [ ( k1Desired(iP(1)) - k1Desired(iP) )' ...
      ( k2Desired(jP(1)) - k2Desired(jP) )' ];
    dist2PeaksPos = sqrt( kVecPos(:,1).^2 + kVecPos(:,2).^2 );
    indsKPos = [ (iP(1) - iP)' (jP(1) - jP)' ];
    ka = mean( dist2PeaksPos( dist2PeaksPos >= 0 ) );
    % real space lattice spacing
    if numPeaksK > 6
      a =  4 * pi / ( sqrt(3) * mean( unique( dist2PeaksPos( dist2PeaksPos >= 0) ) ) ) ;
    elseif numPeaksK == 3
      a = 2 * pi / mean( unique( dist2PeaksPos( dist2PeaksPos >= 0 ) ) );
    else
      a = NaN;
    end
    % plotting
    if plotFlag
      % get real position peaks
      if isfinite( a ) && a < lBox / 2
        % grab some section in the middle
        deltaX1 = ceil( a / lBox * n1 );
        deltaX2 = ceil( a / lBox * n2 );
        if n1 > 1
          inds1 = centerInd1 - deltaX1 :centerInd1 + deltaX1 ;
        else
          inds1 = 1;
        end
        if n2 > 1
          inds2 = centerInd2 - deltaX2 :centerInd2 + deltaX2 ;
        else
          inds2 = 1;
        end
        desiredReal = rho(inds1, inds2);
        newSize = size(desiredReal);
        [~,newMax] = max( desiredReal(:) );
        [iMax, jMax] = ind2sub( newSize , newMax ) ;
        % Center it on a peak
        newCenter1 = inds1(iMax);
        newCenter2 = inds2(jMax);
        if n1 > 1
          inds1 = newCenter1 - deltaX1 : newCenter1 + deltaX1;
        else
          inds1 = 1;
        end
        if n2 > 1
          inds2 = newCenter2 - deltaX2 : newCenter2 + deltaX2;
        else
          inds2 = 1;
        end
        desiredReal = rho(inds1, inds2);
        inds = inds2;
      else
        desiredReal = rho;
        inds = 1:n2;
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
        plot( xVec(inds), desiredReal );
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
  cryPeaks.numPeaksK = numPeaksK;
  cryPeaks.peakValK = peakValK;
  cryPeaks.kVecAll = kVecAll;
  cryPeaks.dist2PeaksAll = dist2PeaksAll;
  cryPeaks.indsKAll = indsKAll;
  cryPeaks.kVecPos = kVecPos;
  cryPeaks.dist2PeaksPos = dist2PeaksPos';
  cryPeaks.indsKPos = indsKPos;
  cryPeaks.dk = dk;
  cryPeaks.a = a;
  cryPeaks.ka = ka;
catch err
  fprintf('%s\n', err.getReport('extended')) ;
  keyboard
end