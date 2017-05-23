function phaseMap = mapPhaseSoftShoulderLinear( bcRange, aRange, rs, lVec, randomVals, numVals, plotPhase )

numLs = length(lVec);
% numVals should be a vector of bc and a
if length( numVals ) == 1
  numVals = [numVals  numVals];
end
% get parameter arrays
if randomVals
  bcVec = linspace( bcRange(1), bcRange(2), max( 100, 10 * numVals(1)  ) );
%   bcVec = bcVec( randperm( length(bcVec), numVals(1) ) );
  aVec = linspace( aRange(1), aRange(2), max( 100, 10 * numVals(2) ) );
%   aVec = aVec( randperm( length(aVec), numVals(2) ) );
  paramMat = combvec( bcVec, aVec);
  numMapParams = size( paramMat, 2 );
  randInds = randperm( numMapParams,  mean(numVals) );
  paramMat = paramMat( :, randInds );
  numMapParams = size( paramMat, 2 );
  bcVec = paramMat(1,:);
  aVec = paramMat(2,:);
  paramMat = combvec( paramMat, lVec );
	numParams = size( paramMat, 2 );
else
  if bcRange(1) == 0
    bcVec = linspace( bcRange(1), bcRange(2), numVals(1)+1 );
    bcVec = bcVec( bcVec > 0 );
  else
    bcVec = linspace( bcRange(1), bcRange(2), numVals(1)+1 );
  end
    aVec = linspace( aRange(1), aRange(2), numVals(2) );
    paramMat = combvec( bcVec, aVec);
    numMapParams = size( paramMat, 2 );
    bcVec = paramMat(1,:);
    aVec = paramMat(2,:);
    paramMat = combvec( bcVec, aVec, lVec );  
end
% build n vec from l
% nVec = ceil( 10 * lVec  / pi + 2 );
% nVec = nVec + mod(nVec,2);
% allocate
phaseMap = zeros( numParams, 1 );
% initalize parameters
for ii = 1:numParams
  bcTemp = paramMat(1,ii);
  aTemp = paramMat(2,ii);
  lTemp = paramMat(3,ii);
  nTemp = ceil( 10 * lTemp  / pi + 2 );
  nTemp = nTemp + mod(nTemp,2);
  paramVec = [ nTemp, nTemp, lTemp, lTemp, 1, aTemp, 1, rs, bcTemp ];
  disp = dispersionSoftShoulder( paramVec );
  phaseMap(ii) = disp.phaseId;
end
% ignore coexistence
phaseMap( phaseMap == 4 ) = 2;
phaseMap( phaseMap == 5 ) = 3;
% reshape to a user-friendly form
phaseMap = reshape( phaseMap, [ numMapParams, numLs ] );
% plot time
if plotPhase
  figure(1000)
  colorArray = colormap(['lines(' num2str(3) ')']);
  close(1000)
  for ii = 1:numLs
    figure()
    tempMat = phaseMap(:,ii);
    scatter( bcVec( tempMat(:) == 1 ), aVec( tempMat(:) == 1 ), ...
      [], colorArray(1,:), 'filled' );
    hold on
    scatter( bcVec( tempMat(:) == 2 ), aVec( tempMat(:) == 2 ), ...
      [], colorArray(2,:), 'filled' );
    scatter( bcVec( tempMat(:) == 3 ), aVec( tempMat(:) == 3 ), ...
      [], colorArray(3,:), 'filled' );
    xlabel('bc'); ylabel('a');
    legend('liquid','crystal A','crystal B', 'location', 'best')
    titstr = ['Linear Phase Diagram L = ' num2str( lVec(ii ) )];
    title(titstr);
  end
end
