% [disp] = dispersionSoftShoulder( paramVec,  plotMe )
%
% get the dispersion relation for soft shoulder potential with
% given input parameter
% paramVec(1) = n1 
% paramVec(2) = n2
% paramVec(3) = l1
% paramVec(4) =l2
% paramVec(5) = lrE1
% paramVec(6) = lrE2
% paramVec(7) = lrL1
% paramVec(8) = lrL2 
% paramVec(9) = c
%
function [disp] = dispersionSoftShoulder( paramVec,  plotMe )
% set plot to zero if no specified
if nargin == 1
  plotMe = 0;
end
% add Subroutine path
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
% move parameters from vec to vals
n1 = paramVec(1);
n2 = paramVec(2);
n3 = 1;
l1 = paramVec(3);
l2 = paramVec(4);
l3 = 1;
lrE1 = paramVec(5);
lrE2 = paramVec(6);
lrL1 = paramVec(7);
lrL2 = paramVec(8);
c = paramVec(9);
% make grid
[gridObj] = GridMakerPBCxk(n1,n2,n3,l1,l2,l3);
% get FT of potential
[~, vFt] = softShoulder2D( lrE1, lrE2, lrL1, lrL2, ...
  n1, l1, n2, l2 );
ck = - l1 * l1 / (n1 .* n2) * real(vFt) ;
% ck = -real(vFt) ;
rho = c;
% dispersion
omega = -( gridObj.k1rep2 .^ 2 + gridObj.k2rep2 .^ 2 ) .* ( 1 - rho .*  ck );
% plot it
center1 = n1/2 + 1;
center2 = n2/2 + 1;
scaleK1 = lrL1 * gridObj.k1;
scaleK2 = lrL1 * gridObj.k2;
maxInd1 = find( scaleK1(center1:end)  > 8 );
maxInd2 = find( scaleK2(center2:end)  > 8 );
inds1 = max(center1-maxInd1,1):min( center1+maxInd1, n1);
inds2 = max(center2-maxInd2,1):min(center2+maxInd2,n1);
indsPos1 = center1:center1+maxInd1;
if isempty(maxInd1) || isempty(maxInd2)
  error( 'Increase the number of gridpoint' )
end
% get slice
% from symmetry, plot a slice ( kx > 0 )
k1pos = scaleK1( indsPos1 );
omegaSlice = omega( indsPos1, center2)';
% get peaks
[omegaPeaks, maxInds] = findpeaks(omegaSlice);
omegaPeaks = [0 omegaPeaks];
maxInds = [1 maxInds];
omegaPeaks = omegaPeaks( omegaPeaks >= 0 );
maxInds = maxInds( omegaPeaks >= 0 );
[omegaMax, omegaMaxInd] = max(omegaPeaks);
kPeaks = k1pos( maxInds );
kPeakMax = kPeaks(omegaMaxInd);
% Guess at the crystal
if length(omegaPeaks) == 1
  phase = 'liquid';
  phaseId = 1;
elseif length(omegaPeaks) == 2
  if kPeakMax < 4
    phase = 'crystal A';
    phaseId = 2;
  else
    phase = 'crystal B';
    phaseId = 3;
  end
elseif length(omegaPeaks) == 3
  phase = '2 unstable peaks: ';
  if kPeakMax < 4
    phase = [phase 'crystal A'];
    phaseId = 4;
  else
    phase = [phase 'crystal B'];
    phaseId = 5;
  end
end
% Store
disp.omega = omega;
disp.omegaSlice = omegaSlice;
disp.omegaPeaks = omegaPeaks;
disp.omegaMax = omegaMax;
disp.kPeaks = kPeaks;
disp.kPeakMax = kPeakMax;
disp.phase = phase;
disp.phaseId = phaseId;
% plots
if plotMe
  % min inds for plotting
  [omegaMin] = findpeaks(-omegaSlice);
  omegaMin = min( -omegaMin );
  figure()
  subplot(1,3,1)
  imagesc( scaleK2(inds2), scaleK1(inds1), ...
    real( rho * ck ( inds1, inds2  )' ) );
  axis square
  xlabel('$$ k_x R $$'); ylabel('$$ k_y R $$');
  title('$$ \rho c(k)$$')
  colorbar
  subplot(1,3,2)
  imagesc( scaleK2, scaleK1, real( omega(inds1, inds2 )' ));
  axis square
  colorbar
  xlabel('$$ k_x R $$'); ylabel('$$ k_y R $$');
  title('$$\omega(k)$$')
  subplot(1,3,3)
  xTemp = k1pos;
  yTemp = omegaSlice;
  p = plot( xTemp, yTemp );
  p.LineWidth = 2;
  hold on
  numPoints = length( k1pos  );
  yTemp = zeros( 1, numPoints );
  p = plot( xTemp, yTemp,'--' );
  p.LineWidth = 2;
  ax = gca;
  ax.XLim = [ 0 max(xTemp) ];
  scalFact  = 0.25; %for limits
  ax.YLim = [ omegaMin + omegaMin * scalFact, max( omegaMax + omegaMax * scalFact, 10 ) ];
  xlabel(' $$ kR $$ ' ); ylabel(' $$ \omega $$ ' );
  axis square
  title('$$\omega(k)$$ slice kx')
end
