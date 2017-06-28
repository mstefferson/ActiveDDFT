% [pDist] = pairDistCalcRho(rho,l1,l2,lRod)
%
% Calculates the density dependent pair correlation function
% g(r,r') = rho^(2) ./ rho^(1)rho^(1)
%
function [pDist] = pairDistCalcRho(rho,l1,l2,lRod,plotflag)
% add paths just in case
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
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
expV =  mayer + 1;
pDist = zeros( n1, n2);
intAng = zeros( n3, n3);
allInds1 = 1:n1;
allInds2 = 1:n2;
allInds3 = 1:n3;
shiftInds1 = -n1/2:n1/2-1;
shiftInds2 = -n2/2:n2/2-1;
% calculate the pair dist
for ii = 1:n1
  shift1Temp = shiftInds1(ii);
  % inds for loop over delta x
  shiftInds1 = mod( (shift1Temp-1) + allInds1 - 1, n1 ) + 1;
  for jj = 1:n2
    shift2Temp = shiftInds2(jj);
    % inds for loop over delta y
    shiftInds2 = mod( (shift2Temp-1) +  allInds2 - 1, n2 ) + 1;
    % integrate of r' each u and u+u'.
    for mm = 1:n3
      shiftInds3 = mod( (mm-1) +  allInds3 - 1, n3 ) + 1;
      for nn = 1:n3
        mat2IntTemp = rho(allInds1, allInds2, allInds3(mm) ) .* rho( shiftInds1, shiftInds2, shiftInds3(nn) );
        tempInt = trapz_periodic( x1, trapz_periodic( x2, mat2IntTemp, 2 ), 1);
        intAng(nn,mm) = expV( ii, jj, shiftInds3(nn), allInds3(mm) ) .* tempInt ;
      end
    end
    pDist( ii, jj ) = trapz_periodic( phi, ...
      trapz_periodic( phi, intAng, 1 ), 2 );
  end
end
% Normalization
nParticles = trapz_periodic( x1, trapz_periodic( x2, trapz_periodic( phi, rho, 3 ), 2 ), 1 );
V = l1 .* l2;
normFac = (nParticles^2 ) ./ V;
pDist = 1 / normFac * pDist;
% center
pDistSave = pDist;
shift1 = round( n1 / 2 );
shift2 = round( n2 / 2 );
pDist = circshift( circshift( pDistSave, shift1, 1 ), shift2, 2 );
% plot
if plotflag
  figure()
  stretchFac = 1.5;
  axisLim = stretchFac * lRod;
  ax1 = subplot(1,2,1);
  imagesc( x1, x2, pDist );
  ax1.XLim = [-axisLim axisLim];
  ax1.YLim = [-axisLim axisLim];
  title('pair distribution g(x,y): new Mayer')
  axis square
  xlabel('x'); ylabel('y')
  ch = colorbar;
  ch.Limits = [0 1];
  ch.Ticks = [0:0.2:1]; 
  % slices
  ax2 = subplot(1,2,2);
  plot( x1, pDist(:,shift2+1), x2, pDist(shift1+1,:) )
  xlabel('position'); ylabel('g(r)');
  title('slice')
  legend('x','y', 'location', 'best' )
  axis square
end
