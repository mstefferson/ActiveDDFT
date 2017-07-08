% [pDist] = pairDistCalcRho(rho,l1,l2,lRod)
%
% Calculates the density dependent pair correlation function
% g(r,r') = rho^(2) ./ rho^(1)rho^(1)
%
function [pDist] = pairDistCalcRho( rho, l1, l2, lRod, plotflag, saveName )
% set saveMe
if nargin == 4
  plotflag = 0;
  saveMe = 0;
elseif nargin == 5
  saveMe = 0;
elseif isempty( saveName )
  saveMe = 0;
else
  saveMe = 1;
end
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
intOverRprime = zeros( n3, n3 );
% set inds
allInds1 = 1:n1;
allInds2 = 1:n2;
allInds3 = 1:n3;
mayerInds1 = [0:n1/2 -n1/2+1:1:-1];
mayerInds2 = [0:n2/2 -n2/2+1:1:-1];
mayerInds3 = [0:n3/2 -n3/2+1:1:-1];
% calculate the pair dist
for ii = 1:n1
  % set x
  r1Temp = mayerInds1(ii); % Actual location of ii 
  % set x+x'
  shiftInds1 = mod( (r1Temp-1) + allInds1 - 1, n1 ) + 1;
  for jj = 1:n2
    % set y
    r2Temp = mayerInds2(jj);
    % set y+y'
    shiftInds2 = mod( (r2Temp-1) +  allInds2 - 1, n2 ) + 1;
    % integrate of r' each u and u+u'.
    for mm = 1:n3
      phi1Temp = mayerInds3(mm);
      shiftInds3 = mod( (phi1Temp-1) +  allInds3 - 1, n3 ) + 1;
      for nn = 1:n3
        mat2IntTemp =  rho( shiftInds1, shiftInds2, shiftInds3(nn) ) .* rho(allInds1, allInds2, allInds3(mm) );
        tempInt = trapz_periodic( x1, trapz_periodic( x2, mat2IntTemp, 2 ), 1);
        intOverRprime(nn,mm) = expV( ii, jj, allInds3(nn), shiftInds3(mm) ) .* tempInt ;
      end
    end
    pDist( ii, jj ) = trapz_periodic( phi, ...
      trapz_periodic( phi, intOverRprime, 1 ), 2 );
  end
  fprintf('%f percent done\n', 100*ii/n1)
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
% Save it
if saveMe 
  save(saveName,'pDist')
end
% plot
if plotflag
  fontSize = 10;
  figure()
  % full
  ax1 = subplot(1,3,1);
  imagesc( x1, x2, pDist' );
  title('pair distribution g(x,y):')
  axis square
  xlabel('x'); ylabel('y')
  ch = colorbar;
  ch.Limits = [0 max(pDist(:))];
  ch.Ticks = 0:0.2:max(pDist(:))+0.2; 
  stretchFac = 1.5;
  axisLim = stretchFac * lRod;
  ax = gca;
  ax.FontSize = fontSize;
  % zoom
  ax2 = subplot(1,3,2);
  imagesc( x1, x2, pDist' );
  ax2.XLim = [-axisLim axisLim];
  ax2.YLim = [-axisLim axisLim];
  title('pair distribution g(x,y):')
  axis square
  xlabel('x'); ylabel('y')
  ch = colorbar;
  ch.Limits = [0 max(pDist(:))];
  ch.Ticks = 0:0.2:max(pDist(:))+0.2; 
  ax = gca;
  ax.FontSize = fontSize;
  % slices
  ax3 = subplot(1,3,3);
  plot( x1, pDist(:,shift2+1), x2, pDist(shift1+1,:) )
  xlabel('position'); ylabel('g(r)');
  title('slice')
  axis square
  leg = legend('x','y');
  leg.Position = [ 0.85 0.45 0.0393 0.05 ];
  leg.FontSize = fontSize -2 ;
  ax = gca;
  ax.FontSize = fontSize - 2;
end
