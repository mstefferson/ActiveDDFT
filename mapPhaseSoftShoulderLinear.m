% phase = mapPhaseSoftShoulderLinear( cRange, aRange, rs, lVec, ...
%  randomVals, numVals, trialId, plotPhase, saveMe )

function phase = mapPhaseSoftShoulderLinear( cRange, aRange, rs, lVec, ...
  randomVals, numVals, trialId, plotPhase, saveMe )
show2Peaks = 1;
% add Subroutine path
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
% grab length of l Vec
numLs = length(lVec);
% get parameter arrays
if randomVals
  cVec = linspace( cRange(1), cRange(2), max( 100, 10 * numVals  ) );
  aVec = linspace( aRange(1), aRange(2), max( 100, 10 * numVals ) );
  paramMat = combvec( cVec, aVec);
  numMapParams = size( paramMat, 2 );
  randInds = randperm( numMapParams,  numVals );
  paramMat = paramMat( :, randInds );
  numMapParams = size( paramMat, 2 );
  cVec = paramMat(1,:);
  aVec = paramMat(2,:);
  paramMat = combvec( paramMat, lVec );
  numParams = size( paramMat, 2 );
else
  % if not random, make a grid of parameter values of length sqrt(numVals)
  numVals = round( sqrt( numVals ) );
  if cRange(1) == 0
    cVec = linspace( cRange(1), cRange(2), numVals+1 );
    cVec = cVec( cVec > 0 );
  else
    cVec = linspace( cRange(1), cRange(2), numVals+1 );
  end
  aVec = linspace( aRange(1), aRange(2), numVals );
  paramMat = combvec( cVec, aVec);
  numMapParams = size( paramMat, 2 );
  cVec = paramMat(1,:);
  aVec = paramMat(2,:);
  paramMat = combvec( cVec, aVec, lVec );
end
% allocate
phaseMat = zeros( numParams, 1 );
% initalize parameters
parfor ii = 1:numParams
  cTemp = paramMat(1,ii);
  aTemp = paramMat(2,ii);
  lTemp = paramMat(3,ii);
  nTemp = ceil( 10 * lTemp  / pi + 2 );
  nTemp = nTemp + mod(nTemp,2);
  paramVec = [ nTemp, nTemp, lTemp, lTemp, 1, aTemp, 1, rs, cTemp ];
  disp = dispersionSoftShoulder( paramVec );
  phaseMat(ii) = disp.phaseId;
end
% reshape to a user-friendly form
phaseMat = reshape( phaseMat, [ numMapParams, numLs ] );
% store it
phase.phaseMap = phaseMat;
phase.cVec = cVec';
phase.aVec = aVec';
phase.rs = rs;
phase.lVec = lVec';
phase.trialId = trialId;
phase.random = randomVals;
phase.numVals = numVals;
% plot time
if plotPhase
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
end
% save it
if saveMe
  filename = ['linstabSoftShoulder' ...
    '_c' num2str( cRange(1) ) '-' num2str( cRange(2) ) ...
    '_a' num2str( aRange(1) ) '-' num2str( aRange(2) ) ...
    '_lnum' num2str( length(l) )  ...
    '_tr' num2str( trialId ), '.mat' ];
    save( filename, 'phase');
    if ~exist('linstab', 'dir'); mkdir('./linstab'); end
  movefile( filename, './linstab');
end

