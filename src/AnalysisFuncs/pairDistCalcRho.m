% [pDist] = pairDistCalcRho( rho, l1, l2, lRod, plotflag, saveName )
%
% Calculates the density dependent pair correlation function
% g(r,r') = rho^(2) ./ rho^(1)rho^(1)
%
function [pDist] = pairDistCalcRho( rho, l1, l2, lRod, ...
  indsWant1, indsWant2, calcG1G2Flag, ...
  plotflag, saveName )
% get sizes
[n1,n2,n3] = size( rho );
% set saveMe
if nargin == 4
  calcG1G2Flag = 0;
  indsWant1 = 1:n1;
  indsWant2 = 1:n2;
  plotflag = 0;
  saveMe = 0;
elseif nargin == 5
  indsWant2 = 1:n2;
  calcG1G2Flag = 0;
  plotflag = 0;
  saveMe = 0;
elseif nargin == 6
  calcG1G2Flag = 0;
  plotflag = 0;
  saveMe = 0;
elseif nargin == 7
  plotflag = 0;
  saveMe = 0;
elseif isempty( saveName )
  plotflag = 0;
  saveMe = 0;
else
  saveMe = 1;
end
% get inds to integrate
inds1 = intersect( indsWant1, 1:n1 );
inds2 = intersect( indsWant2, 1:n2 );

%[inds1, inds2] =  combvec( inds1, inds2 );
numInds1 = length( inds1 );
numInds2 = length( inds2 );

% build grid
dx1 = l1 ./ n1;
dx2 = l2 ./ n2;
dphi = 2 * pi / n3;
x1 = dx1 .* (-n1/2:n1/2-1 );
x2 = dx2 .* (-n2/2:1:n2/2-1 );
phi = 0 : dphi : 2*pi - dphi;
% new mayer: integrate over angles from the start
[mayer] = mayerFncHrLabFrame( n1, n2, n3, l1, l2, lRod );
% initialize intergrator class
pIntegrator = PairDistIntegratorClass( n1, n2, n3, calcG1G2Flag, x1, x2, phi, mayer );
delta0 = zeros( numInds1, numInds2); % for pair dist
delta1 = zeros( numInds1, numInds2); % for polar dist
delta2 = zeros( numInds1, numInds2); % for nem dist
% track progress
ticId = tic;
trackProgMod =  ceil( n1 / 100);
% calculate the pair dist
for ii = 1:numInds1
  % update shifted inds
  pIntegrator.updateShiftInds( inds1(ii),  1 );
  for jj = 1:numInds2
    % update shifted inds
    pIntegrator.updateShiftInds( inds2(jj),  2 );
    pIntegrator.calcDeltaIntegrals(rho);
    delta0( ii, jj ) = pIntegrator.Delta0;
    if calcG1G2Flag
      delta1( ii, jj ) = pIntegrator.Delta1;
      delta2( ii, jj ) = pIntegrator.Delta2;
    end
  end
  if mod( ii, trackProgMod  ) == 0
    fprintf('%f percent done\n', 100*ii/numInds1)
  end
end
timeRun = toc(ticId);
fprintf('runtime: %.1f sec\n', timeRun);
% Normalization
nParticles = trapz_periodic( x1, trapz_periodic( x2, trapz_periodic( phi, rho, 3 ), 2 ), 1 );
V = l1 .* l2;
normFac = (nParticles^2 ) ./ V;
% calculate pair distributions
pDist0 = 1 / normFac * delta0; % pair dist
if calcG1G2Flag
  pDist1 = delta1 ./ delta0; % polar dist
  pDist1(1,1) = 0; % get rid of dividing by zero
  pDist2 = (2*delta2-delta0) ./ delta0; % nem dist
  pDist2(1,1) = 0; % get rid of dividing by zero
end
% Rotate and center it
shiftColumn = round( n1 / 2 );
shiftRow = round( n2 / 2 );
pDist0Center = circshift( pDist0, [round( n1 / 2 ) round( n2 / 2 ) ] );
% center it. rows and columns are flopped from rotating
if calcG1G2Flag
  pDist1Center = circshift( pDist1, [round( n1 / 2 ) round( n2 / 2 ) ] );
  pDist2Center = circshift( pDist2, [round( n1 / 2 ) round( n2 / 2 ) ] );pDist2RotCenter = circshift( circshift( pDist2Rot, -shiftRow+1, 1 ), shiftColumn, 2 );
end
% store it
pDist.pDist0 = pDist0;
pDist.pDist0Center = pDist0Center;
if calcG1G2Flag
  pDist.pDist1 = pDist1;
  pDist.pDist2 = pDist2;
  pDist.pDist1Center = pDist1Center;
  pDist.pDist2Center = pDist2Center;
end
% Save it
if saveMe
  save(saveName,'pDist')
end

% plot
if plotflag
  % plot p0
  ttlstr = '$$g_0(x,y)$$';
  pairDistPlotSingle( pDist0RotCenter, l2, l1, ttlstr);
  % plot the three order parameter distribution functions
  if calcG1G2Flag
    pairDistPlotOps( pDist0RotCenter, pDist1RotCenter, pDist2RotCenter, l2, l1)
  end
end
