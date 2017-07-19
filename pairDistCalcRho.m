% [pDist] = pairDistCalcRho( rho, l1, l2, lRod, plotflag, saveName )
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
delta0 = zeros( n1, n2); % for pair dist
delta1 = zeros( n1, n2); % for polar dist
delta2 = zeros( n1, n2); % for nem dist
intOverRprime0= zeros( n3, n3 );
intOverRprime1= zeros( n3, n3 );
intOverRprime2= zeros( n3, n3 );
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
      allInds3mm = allInds3(mm);
      for nn = 1:n3
        shiftInds3nn = shiftInds3(nn);
        mat2Intdelta0 =  rho( shiftInds1, shiftInds2, shiftInds3nn ) .* rho(allInds1, allInds2, allInds3mm );
        mat2Intdelta1 =  cos(phi1Temp) * rho( shiftInds1, shiftInds2, shiftInds3nn ) .* rho(allInds1, allInds2, allInds3mm );
        mat2Intdelta2 =  cos(phi1Temp)^2 * rho( shiftInds1, shiftInds2, shiftInds3nn ) .* rho(allInds1, allInds2, allInds3mm );
        tempInt0 = trapz_periodic( x1, trapz_periodic( x2, mat2Intdelta0, 2 ), 1);
        tempInt1 = trapz_periodic( x1, trapz_periodic( x2, mat2Intdelta1, 2 ), 1);
        tempInt2 = trapz_periodic( x1, trapz_periodic( x2, mat2Intdelta2, 2 ), 1);
        intOverRprime0(nn,mm) = expV( ii, jj, allInds3(nn), shiftInds3(mm) ) .* tempInt0 ;
        intOverRprime1(nn,mm) = expV( ii, jj, allInds3(nn), shiftInds3(mm) ) .* tempInt1 ;
        intOverRprime2(nn,mm) = expV( ii, jj, allInds3(nn), shiftInds3(mm) ) .* tempInt2 ;
      end
    end
    delta0( ii, jj ) = trapz_periodic( phi, ...
      trapz_periodic( phi, intOverRprime0, 1 ), 2 );
    delta1( ii, jj ) = trapz_periodic( phi, ...
      trapz_periodic( phi, intOverRprime1, 1 ), 2 );
    delta2( ii, jj ) = trapz_periodic( phi, ...
      trapz_periodic( phi, intOverRprime2, 1 ), 2 );
  end
  fprintf('%f percent done\n', 100*ii/n1)
end
% Normalization
nParticles = trapz_periodic( x1, trapz_periodic( x2, trapz_periodic( phi, rho, 3 ), 2 ), 1 );
V = l1 .* l2;
normFac = (nParticles^2 ) ./ V;
% calculate pair distributions
pDist0 = 1 / normFac * delta0; % pair dist
pDist1 = delta1 ./ delta0; % polar dist
pDist1(1,1) = 0; % get rid of dividing by zero
pDist2 = (2*delta2-delta0) ./ delta0; % nem dist
pDist2(1,1) = 0; % get rid of dividing by zero
% Rotate it 
pDist0Rot = rot90(pDist0);
pDist1Rot = rot90(pDist1);
pDist2Rot = rot90(pDist2);
% center it. rows and columns are flopped from rotating 
shiftColumn = round( n1 / 2 ); 
shiftRow = round( n2 / 2 );
pDist0RotCenter = circshift( circshift( pDist0Rot, -shiftRow+1, 1 ), shiftColumn, 2 );
pDist1RotCenter = circshift( circshift( pDist1Rot, -shiftRow+1, 1 ), shiftColumn, 2 );
pDist2RotCenter = circshift( circshift( pDist2Rot, -shiftRow+1, 1 ), shiftColumn, 2 );
% store it
pDist.pDist0 = pDist0;
pDist.pDist1 = pDist1;
pDist.pDist2 = pDist2;
pDist.pDist0RotCenter = pDist0RotCenter;
pDist.pDist1RotCenter = pDist1RotCenter;
pDist.pDist2RotCenter = pDist2RotCenter;
% Save it
if saveMe 
  save(saveName,'pDist')
end

% plot
if plotflag
  % plot p0
  ttlstr = '$$g_0(x,y)$$';
  pairDistPlot( pDist0RotCenter, l2, l1, ttlstr)
  % plot the three pair distributions
  figure()
  % g0
  maxVal = max(pDist0(:));
  minVal = min( min(pDist0(:)), 0 );
  ax0 = subplot(1,3,1);
  ax0.FontSize = fontSize;
  imagesc( x1, x2, pDist0RotCenter );
  title('$$g_0(x,y)$$')
  axis square
  xlabel('$$x$$'); ylabel('$$y$$')
  ax0.XLim = [-xLim xLim];
  ax0.YLim = [-yLim yLim];
  ax0.CLim = [minVal maxVal];
  ax0.YDir = 'normal';
  colorbar;
  % g1 (polar)
  ax1 = subplot(1,3,2);
  maxVal = max(pDist1(:));
  minVal = min( min(pDist1(:)), 0 );
  ax1.FontSize = fontSize;
  imagesc( x1, x2, pDist1RotCenter );
  title('$$ g_1(x,y) $$')
  axis square
  xlabel('$$x$$'); ylabel('$$y$$')
  ax1.XLim = [-xLim xLim];
  ax1.YLim = [-yLim yLim];
  ax1.CLim = [minVal maxVal];
  ax1.YDir = 'normal';
  colorbar;
  % g2 (nem)
  maxVal = max(pDist2(:));
  minVal = min( min(pDist2(:)), 0 );
  ax2 = subplot(1,3,3);
  ax2.FontSize = fontSize;
  imagesc( x1, x2, pDist2RotCenter );
  title('$$ g_2 (x,y) $$')
  axis square
  xlabel('$$x$$'); ylabel('$$y$$')
  ax2.XLim = [-xLim xLim];
  ax2.YLim = [-yLim yLim];
  ax2.CLim = [minVal maxVal];
  ax2.YDir = 'normal';
  colorbar;
end

