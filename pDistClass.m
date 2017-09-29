% [pDist] = pairDistCalcRho( rho, l1, l2, lRod, plotflag, saveName )
%
% Calculates the density dependent pair correlation function
% g(r,r') = rho^(2) ./ rho^(1)rho^(1)
%
function [pDist] = pDistClass( rho, l1, l2, lRod, calcG1G2Flag, ...
  plotflag, saveName )
% set saveMe
if nargin == 4
  calcG1G2Flag = 0;
  plotflag = 0;
  saveMe = 0;
elseif nargin == 5
  plotflag = 0;
  saveMe = 0;
elseif nargin == 6
  plotflag = 0;
  saveMe = 0;
elseif isempty( saveName )
  saveMe = 0;
else
  saveMe = 1;
end
% build phi
[n1,n2,n3] = size( rho );
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
delta0 = zeros( n1, n2); % for pair dist
delta1 = zeros( n1, n2); % for polar dist
delta2 = zeros( n1, n2); % for nem dist
% track progress
ticId = tic;
trackProgMod =  ceil( n1 / 100);
% calculate the pair dist
for ii = 1:n1
  % update shifted inds
  pIntegrator.updateShiftInds( ii,  1 );
  for jj = 1:n2
    % update shifted inds
    pIntegrator.updateShiftInds( jj,  2 );
    pIntegrator.calcDeltaIntegrals(rho);
    delta0( ii, jj ) = pIntegrator.Delta0;
    if calcG1G2Flag
      delta1( ii, jj ) = pIntegrator.Delta1;
      delta2( ii, jj ) = pIntegrator.Delta2;
    end
  end
  if mod( ii, trackProgMod  ) == 0
    fprintf('%f percent done\n', 100*ii/n1)
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
pDist0Rot = rot90(pDist0);
pDist0RotCenter = circshift( circshift( pDist0Rot, -shiftRow+1, 1 ), shiftColumn, 2 );
% center it. rows and columns are flopped from rotating
if calcG1G2Flag
  pDist1Rot = rot90(pDist1);
  pDist2Rot = rot90(pDist2);
  pDist1RotCenter = circshift( circshift( pDist1Rot, -shiftRow+1, 1 ), shiftColumn, 2 );
  pDist2RotCenter = circshift( circshift( pDist2Rot, -shiftRow+1, 1 ), shiftColumn, 2 );
end
% store it
pDist.pDist0 = pDist0;
pDist.pDist0RotCenter = pDist0RotCenter;
if calcG1G2Flag
  pDist.pDist1 = pDist1;
  pDist.pDist2 = pDist2;
  pDist.pDist1RotCenter = pDist1RotCenter;
  pDist.pDist2RotCenter = pDist2RotCenter;
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
