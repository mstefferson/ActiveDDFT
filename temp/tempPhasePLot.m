%%
show2Peaks = 1;
phaseMat = phase.phaseMap;
aRange = phase.aRange;
cRange = phase.cRange;
numLs = length(phase.lVec);
cVec = phase.cVec;
aVec = phase.aVec;
lVec = phase.lVec;
%%
% ignore coexistence
if show2Peaks
  legcell = {'liquid', 'crystal A', 'crystal B', ...
    '2 peaks: A', '2 peaks: B'};
  numColors = 5;
else
  phaseMat( phaseMat == 4 ) = 2;
  phaseMat( phaseMat == 5 ) = 3;
  legcell = {'liquid', 'crystal A', 'crystal B'};
  numColors = 3;
end
figure(1000)
colorArray = colormap(['lines(' num2str(numColors) ')']);
close(1000)
for ii = 1:numLs
  figure()
  tempMat = phaseMat(:,ii);
  scatter( cVec( tempMat(:) == 1 ), aVec( tempMat(:) == 1 ), ...
    [], colorArray(1,:), 'filled' );
  hold on
  scatter( cVec( tempMat(:) == 2 ), aVec( tempMat(:) == 2 ), ...
    [], colorArray(2,:), 'filled' );
  scatter( cVec( tempMat(:) == 3 ), aVec( tempMat(:) == 3 ), ...
    [], colorArray(3,:), 'filled' );
  if show2Peaks
    scatter( cVec( tempMat(:) == 4 ), aVec( tempMat(:) == 4 ), ...
      [], colorArray(4,:), 'filled' );
    scatter( cVec( tempMat(:) == 5 ), aVec( tempMat(:) == 5 ), ...
      [], colorArray(5,:), 'filled' );
  end
  % clean up plot
  axis square
  xlabel('$$\rho R^2 $$'); ylabel('$$ a $$');
  leg = legend( legcell );
  leg.Position = [0.7184 0.4672 0.1098 0.2118];
  ax = gca;
  ax.XTick = cRange(1):(cRange(2)-cRange(1))/4:cRange(2);
  ax.YTick = aRange(1):(aRange(2)-aRange(1))/4:aRange(2);
  ax.XLim = [cRange(1) cRange(2)];
  ax.YLim = [aRange(1) aRange(2)];
  titstr = ['Linear Phase Diagram L = ' num2str( lVec(ii ) )];
  title(titstr);
end
