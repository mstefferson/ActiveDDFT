% [disp] = dispersionSoftShoulder( paramVec,  plotMe )
%
% get the dispersion relation for soft shoulder potential with
% given input parameter
% paramVec(1) = n1
% paramVec(2) = n2
% paramVec(3) = n3
% paramVec(4) = l1
% paramVec(5) = l2
% paramVec(6) = l3
% paramVec(7) = lrE1
% paramVec(8) = lrE2
% paramVec(9) = lrL1
% paramVec(10) = lrL2
% paramVec(11) = c
% paramVec(12) = vD
%
function [stabInfo] = stabilityDrivenSoftIso( paramVec, plotMe )
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
n3 = paramVec(3);
l1 = paramVec(4);
l2 = paramVec(5);
l3 = paramVec(6);
lrE1 = paramVec(7);
lrE2 = paramVec(8);
lrL1 = paramVec(9);
lrL2 = paramVec(10);
c = paramVec(11);
vD = paramVec(12);
% make grid
[gridObj] = GridMakerPBCxk(n1,n2,n3,l1,l2,l3);
% get FT of potential
[~, vFt] = softShoulder2D( lrE1, lrE2, lrL1, lrL2, ...
  n1, l1, n2, l2 );
ck = - l1 * l1 / (n1 .* n2) * real(vFt) ;
% loop over all k
k3 = gridObj.k3;
omegaRealFull = zeros( n1 , n2);
omegaImagFull = zeros( n1 , n2);
deltaM0 = zeros(1,n3);
deltaM0(n3/2+1) = 1;
% Full thing
for ii = 1:n1
  for jj = 1:n2
    k1 = gridObj.k1(ii);
    k2 = gridObj.k2(jj);
    ckFixedK = ck(ii,jj) * deltaM0;
    omegaDiag = -( k1 .^ 2 + k2 .^ 2 + k3 .^ 2) .* ( 1 - c .*  ckFixedK );
    p1Prefac = -vD / 2 * ( sqrt(-1) * k1 - k2 );
    m1Prefac = -vD / 2 * ( sqrt(-1) * k1 + k2 );
    omegaP1 = p1Prefac * ones(n3,1);
    omegaM1 = m1Prefac * ones(n3,1);
    % build the mat
    omegaMat = diag( omegaDiag' ) + ...
      diag(  omegaP1(1:n3-1) ,1) + diag(  omegaM1(1:n3-1) ,-1);
    % periodic boundaries
    omegaMat(1,n3) = p1Prefac;
    omegaMat(n3,1) = m1Prefac;
    [eigVals] = eig( omegaMat );
    [maxEig, maxInds] = max( real(eigVals) );
    omegaImagFull(ii,jj) = imag( eigVals(maxInds) );
    omegaRealFull(ii,jj) = maxEig ;
  end
end
% store it
stabInfo.omegaRealFull = omegaRealFull;
stabInfo.omegaImagFull = omegaImagFull;
stabInfo.k1 = gridObj.k1;
stabInfo.k2 = gridObj.k2;
% plot if you want
if plotMe
  % full dispersion
  centerk1 = n1/2+1;
  centerk2 = n2/2+1;
  figure()
  subplot(1,3,1)
  imagesc( k2, k1, omegaRealFull' )
  ax = gca;
  ax.YDir = 'normal';
  axis square
  colorbar
  xlabel('$$ k_x R $$'); ylabel('$$ k_y R $$');
  title('$$\max( Re( \omega ) )$$')
  subplot(1,3,2)
  imagesc( k2, k1, omegaImagFull' )
  ax = gca;
  ax.YDir = 'normal';
  axis square
  colorbar
  xlabel('$$ k_x R $$'); ylabel('$$ k_y R $$');
  title('$$\max( Im( \omega ) )$$')
  subplot(1,3,3)
  % compare vs no driving
  paramVec = [n1 n2 l1 l2 lrE1 lrE2 lrL1 lrL2 c];
  [disper] = dispersionSoftShoulder( paramVec,  0);
  inds = centerk1:n1;
  plot( gridObj.k1(inds), omegaRealFull(inds,centerk2), ...
    gridObj.k1(inds), disper.omega(inds,centerk2) )
   axis square
     xlabel('$$ k_x R $$'); ylabel('$$\max( Re( \omega ) )$$');
  legend( 'driving', 'no driving')
end
