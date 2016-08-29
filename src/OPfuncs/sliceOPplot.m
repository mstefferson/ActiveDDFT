% Test function to plot slices of the order parameters

function sliceOPplot( cFinal, poFinal, noFinal, systemObj, ...
  gridObj, rhoFinal, saveTag )

ampFract  = 0.05; % for pos, phi plot
Nx = systemObj.Nx;
Ny = systemObj.Ny;
Nm = systemObj.Nm;
x  = gridObj.x;
y  = gridObj.y;
phi = gridObj.phi;

rowVar = std( cFinal( :, 1 ) );
colVar = std( cFinal( 1, : ) );

% Take slices
if rowVar >= colVar
  rows = 1:Nx;
  cols = 1;
  height = 1:Nm;
  NposVar = Nx;
  posVar = x;
  shiftDim = 1;
  centerPos = Nx / 2;
  
else
  rows = 1;
  cols = 1:Ny;
  height = 1:Nm;
  NposVar = Ny;
  posVar = y;
  posVarLab = 'y';
  shiftDim = 2;
  centerPos = Ny / 2;
end

[maxC, maxCind] = max( cFinal( rows,cols ) );
[maxPO, maxPOind] = max( poFinal( rows,cols ) );
[maxNO, ~] = max( noFinal( rows,cols ) );

% Shift things to the center
cFinal = circshift( cFinal, centerPos - maxCind, shiftDim  );
poFinal = circshift( poFinal, centerPos - maxCind, shiftDim  );
noFinal = circshift( noFinal, centerPos - maxCind, shiftDim  );

% distro slice
distroSlice = ...
  reshape( rhoFinal( rows, cols, height ), [Nx Nm] );
% Shift it
distroSlice = circshift( distroSlice, centerPos - maxCind, 1);
maxDist = max(max( distroSlice ) );
distroSliceChop = distroSlice;
distroSliceChop( distroSlice >  maxDist * ampFract ) = 0;

% Plot various distros
deltaCPpeak = abs( maxCind - maxPOind );
deltaInd = ...
  [-2*deltaCPpeak  -deltaCPpeak 0 deltaCPpeak 2*deltaCPpeak];
plotInd = mod( centerPos + deltaInd -1 , Nm ) + 1;

% Plot OPs
figure()
subplot(1,3,1)
imagesc( x, y, cFinal' );
axis square
colorbar
title( 'C' );
xlabel('x'); ylabel('y');
Ax = gca;
Ax.YDir = 'normal';

subplot(1,3,2)
imagesc( x, y, poFinal' );
axis square
colorbar
title( 'PO' );
xlabel('x'); ylabel('y');
Ax = gca;
Ax.YDir = 'normal';
Ax.CLim = [0 1];

subplot(1,3,3)
imagesc( x, y, noFinal' );
axis square
colorbar
title( 'NO' );
xlabel('x'); ylabel('y');
Ax = gca;
Ax.YDir = 'normal';
Ax.CLim = [0 1];

% Save it
tempTag = '_OPsurf';
savefig( gcf, [saveTag tempTag '.fig'] );
saveas( gcf, [saveTag tempTag '.jpg'], 'jpg')

% Plot a slide through space
figure()
subplot(2,2,1)
plot( posVar, cFinal( rows, cols )./ maxC,...
  posVar, poFinal( rows, cols )./ maxPO, ...
  posVar, noFinal( rows, cols )./ maxNO );
xlabel(posVarLab); ylabel('Normalized amp');
title('Normalized Order Parameters')
legend('c','po', 'no' )

subplot(2,2,2)
plot( posVar, cFinal( rows, cols ) );
xlabel(posVarLab); ylabel('C')
title('C')
Ax = gca;
Ax.YDir = 'normal';

subplot(2,2,3)
plot( posVar, poFinal( rows, cols ) );
xlabel(posVarLab); ylabel('PO')
title('PO')
Ax = gca;
Ax.YDir = 'normal';

subplot(2,2,4)
plot( posVar, noFinal( rows, cols ) );
xlabel(posVarLab); ylabel('NO')
title('NO')
Ax = gca;
Ax.YDir = 'normal';

% Save it
tempTag = '_OPslice';
savefig( gcf, [saveTag tempTag '.fig'] );
saveas( gcf, [saveTag tempTag '.jpg'], 'jpg')

% Plot distribution at different positions
figure()
subplot(1,2,1)
pcolor(posVar, phi,distroSlice)
axis square
colorbar
xlabel('\phi'); ylabel(posVarLab);
title('Distibution sliced through position and angle')
Ax = gca;
Ax.YDir = 'normal';
shading interp;

subplot(1,2,2)
pcolor(posVar, phi,distroSliceChop)
axis square
colorbar
xlabel('\phi'); ylabel(posVarLab);
title('Distribution with high ampiltude sliced out')
Ax = gca;
Ax.YDir = 'normal';
shading interp;

% Save it
tempTag = '_DistroSurf';
savefig( gcf, [saveTag tempTag '.fig'] );
saveas( gcf, [saveTag tempTag '.jpg'], 'jpg')

figure()
counter = 1;
for ii = 1:length( plotInd )  
  if plotInd(ii) ~= centerPos
    subplot( 2,2, counter )
    plot( phi,  distroSlice( plotInd(ii), : ), ...
      phi,  distroSlice( centerPos, : ))
    counter = counter + 1;
    titStr = sprintf('Distro at fixed pos (Cmax - POmax) = %.1f', ...
      deltaCPpeak ./ NposVar );
    title(titStr)
    xlabel('\phi'); ylabel('f(\phi)')
    
    deltaStr = sprintf( 'Cmax + (%.1f) ', deltaInd(ii) ./ NposVar );
    legend( deltaStr, 'C max' );
  end  
end

% Save it
tempTag = '_DistroSlices';
savefig( gcf, [saveTag tempTag '.fig'] );
saveas( gcf, [saveTag tempTag '.jpg'], 'jpg')

