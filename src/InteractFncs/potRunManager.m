function [potObj] =  potRunManager( pot )
% Build potObj
numPot = length(pot);
if numPot
  % Allocate
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
  [ numRuns, ~] = size(paramMat');
  potInds = 1:numRuns;
  paramSepCell = mat2cell( paramMat', ones(1, numRuns), potInputs );
  % Reallocate
  runTemp = cell( numPot, 1);
  potParams = cell( numRuns, 1);
  potStr = cell( numRuns, 1);
  for ii = 1:numRuns
    strTemp = ['_'];
    for jj = 1:numPot
      runTemp{jj} = { names{jj} paramSepCell{ii,jj} };
      strTemp = [  strTemp names{jj} num2str( paramSepCell{ii,jj}, '_%.2f') '_' ];
    end
    potParams{ii} = runTemp;
    potStr{ii} = strTemp;
  end 
else
  potParams = {[]};
  potStr = {''};
  potInds = 1;
end
% Store it
potObj.param = potParams;
potObj.str = potStr;
potObj.inds = potInds;
