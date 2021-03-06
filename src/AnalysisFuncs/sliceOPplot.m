% Test function to plot slices of the order parameters

function sliceOPplot( cFinal, poFinal, noFinal, systemObj, ...
  gridObj, rhoFinal, saveTag )

if nargin == 6
  saveFlag = 0;
else
  saveFlag = 1;
end

varCutoff = 1e-6; % If variance is below this, edits axis;
axisFactor =  1e3;
ampFract  = 0.05; % for pos, phi plot
n1 = systemObj.n1;
n2 = systemObj.n2;
n3 = systemObj.n3;
x  = gridObj.x1;
y  = gridObj.x2;
phi = gridObj.x3;

rowVar = std( cFinal( :, 1 ) );
colVar = std( cFinal( 1, : ) );

% Take slices
if rowVar >= colVar
  maxVar = rowVar;
  rows = 1:n1;
  cols = 1;
  height = 1:n3;
  NposVar = n1;
  posVar = x;
  posVarLab = 'x';
  shiftDim = 1;
  centerPos = n1 / 2;
else
  maxVar = colVar;
  rows = 1;
  cols = 1:n2;
  height = 1:n3;
  NposVar = n2;
  posVar = y;
  posVarLab = 'y';
  shiftDim = 2;
  centerPos = n2 / 2;
end

if maxVar < varCutoff
  aveC    = mean( mean( cFinal ) );
  deltaC  = aveC / axisFactor;
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
  reshape( rhoFinal( rows, cols, height ), [NposVar n3] );
% Shift it
distroSlice = circshift( distroSlice, centerPos - maxCind, 1);
maxDist = max(max( distroSlice ) );
distroSliceChop = distroSlice;
distroSliceChop( distroSlice >  maxDist * ampFract ) = 0;
if maxVar < varCutoff
  aveDistro    = mean( mean( distroSlice ) );
  deltaDistro  = aveDistro / axisFactor;
end

% Plot various tros
deltaCPpeak = abs( maxCind - maxPOind );
deltaInd = ...
  [-2*deltaCPpeak  -deltaCPpeak 0 deltaCPpeak 2*deltaCPpeak];
plotInd = mod( centerPos + deltaInd -1 , NposVar ) + 1;

% Plot OPs
figure()
subplot(1,3,1)
imagesc( x, y, cFinal' );
axis square
colorbar
title( 'C' );
xlabel('x'); ylabel('y');
Ax = gca;
Ax.YDir = 'rev';

if maxVar < varCutoff
  if aveC == 0
   Ax.CLim = [ 0 0.1]; 
  else
  Ax.CLim = [ aveC-deltaC aveC+deltaC]; 
  end
end

subplot(1,3,2)
imagesc( x, y, poFinal' );
axis square
colorbar
title( 'PO' );
xlabel('x'); ylabel('y');
Ax = gca;
Ax.YDir = 'rev';
Ax.CLim = [0 1];

subplot(1,3,3)
imagesc( x, y, noFinal' );
axis square
colorbar
title( 'NO' );
xlabel('x'); ylabel('y');
Ax = gca;
Ax.YDir = 'rev';
Ax.CLim = [0 1];

% Save it
if saveFlag
  tempTag = '_OPsurf';
  savefig( gcf, [saveTag tempTag '.fig'] );
  saveas( gcf, [saveTag tempTag '.jpg'], 'jpg')
end

% Plot a slide through space
figure()
subplot(2,2,1)
plot( posVar, cFinal( rows, cols )./ maxC,...
  posVar, poFinal( rows, cols )./ maxPO, ...
  posVar, noFinal( rows, cols )./ maxNO );
xlabel(posVarLab); ylabel('Normalized amp');
title('Normalized Order Parameters')
legend('C','P', 'N' )

subplot(2,2,2)
plot( posVar, cFinal( rows, cols ) );
xlabel(posVarLab); ylabel('C')
title('C')
Ax = gca;
Ax.YDir = 'normal';
if maxVar < varCutoff
  if aveC == 0
   Ax.CLim = [ 0 0.1]; 
  else
  Ax.CLim = [ aveC-deltaC aveC+deltaC]; 
  end
end

subplot(2,2,3)
plot( posVar, poFinal( rows, cols ) );
xlabel(posVarLab); ylabel('P')
title('Polar Order')
Ax = gca;
Ax.YDir = 'normal';

subplot(2,2,4)
plot( posVar, noFinal( rows, cols ) );
xlabel(posVarLab); ylabel('N')
title('Nematic Order')
Ax = gca;
Ax.YDir = 'normal';

% Save it
if saveFlag
  tempTag = '_OPslice';
  savefig( gcf, [saveTag tempTag '.fig'] );
  saveas( gcf, [saveTag tempTag '.jpg'], 'jpg')
end

% Plot distribution at different positions
figure()
subplot(1,2,1)
imagesc(posVar, phi,distroSlice)
axis square
colorbar
xlabel(posVarLab); ylabel('$$\phi$$');
title('Distibution sliced through position and angle')
Ax = gca;
Ax.YDir = 'normal';
shading interp;
if maxVar < varCutoff
  if aveC == 0
   Ax.CLim = [ 0 0.1]; 
  else
  Ax.CLim =  [ aveDistro-deltaDistro aveDistro+deltaDistro]; 
  end
end
subplot(1,2,2)
imagesc(posVar, phi,distroSliceChop)
axis square
colorbar
xlabel(posVarLab); ylabel('$$\phi$$');
title('Distribution with high ampiltude sliced out')
Ax = gca;
Ax.YDir = 'normal';
shading interp;

if maxVar < varCutoff
  if aveC == 0
   Ax.CLim = [ 0 0.1]; 
  else
  Ax.CLim =  [ aveDistro-deltaDistro aveDistro+deltaDistro]; 
  end
end
% Save it
if saveFlag
  tempTag = '_DistroSurf';
  savefig( gcf, [saveTag tempTag '.fig'] );
  saveas( gcf, [saveTag tempTag '.jpg'], 'jpg')
end

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
    xlabel('$$\phi$$'); ylabel('$$f(\phi)$$')
    if maxVar < varCutoff; ...
        Ax.YLim = [ aveDistro-deltaDistro aveDistro+deltaDistro]; 
    end
    deltaStr = sprintf( 'Cmax + (%.1f) ', deltaInd(ii) ./ NposVar );
    legend( deltaStr, 'C max' );
  end
end
% Save it
if saveFlag
  tempTag = '_DistroSlices';
  savefig( gcf, [saveTag tempTag '.fig'] );
  saveas( gcf, [saveTag tempTag '.jpg'], 'jpg')
end

