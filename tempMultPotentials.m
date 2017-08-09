pot = { {'pot1', [1.1, 2.1 3.1], 4.1}; ...
  {'pot2', 1.2, [2.2 3.2], 4.2};...
  {'pot3', [5.3 6.3] } };
% pot = {  };
if ~isempty(pot)
  numPot = length(pot);
  potInputs = zeros( 1, numPot );
  names = cell( 1, numPot );
  for ii = 1:numPot
    potTemp = pot{ii};
    potInputs(ii) = length(potTemp) - 1;
    paramSingMat = potTemp{2};
    names{ii} = potTemp{1};
    for jj = 3:potInputs(ii)+1
      paramSingMat = combvec(paramSingMat,potTemp{jj});
    end
    if ii == 1
      paramMat = paramSingMat;
    else
      paramMat = combvec(paramMat, paramSingMat);
    end
  end
  % put transpose in cell
  [ numRuns, paramsPerRun] = size(paramMat');
  potInds = 1:length(numRuns);
  paramSepCell = mat2cell( paramMat', ones(1, numRuns), potInputs );
  % string
  strCell = cell( numRuns, numPot );
  temp = cell2mat( paramSepCell(:,1) ) ;
  runTemp = cell( numPot, 1);
  runObjCell = cell( numRuns, 1);
  runStrCell = cell( numRuns, 1);
  for ii = 1:numRuns
    strTemp = ['_'];
    for jj = 1:numPot
      runTemp{jj} = { names{jj} paramSepCell{ii,jj} };
      strTemp = [  strTemp names{jj} num2str( paramSepCell{ii,jj}, '_%.2f') '_' ];
    end
    runObjCell{ii} = runTemp;
    runStrCell{ii} = strTemp;
  end 
else
  runObjCell = {[]};
  runStrCell = {''};
  potInds = 1;
end





