% long range interaction and external potential manager
function [potObj] =  potRunManager( pot, interactPotFlag )
% Build potObj
numPot = length(pot);
if numPot
  % number of strings
  % Allocate
  potInputs = zeros( 1, numPot );
  names = cell( 1, numPot );
  correlation = cell( 1, numPot );
  for ii = 1:numPot
    potTemp = pot{ii};
    names{ii} = potTemp{1};
    % if no correlation function given, assume MF
    if interactPotFlag
      if isstr(potTemp{2})
        if ~strcmp( potTemp{2}, 'mf' ) || ~strcmp( potTemp{2}, 'mf' ) 
          fprintf('Erroneous correlation, setting to mf \n')
           potTemp{2} = 'mf';
           correlation{ii} = 'mf';
        end
        strNum = 2;
        paramInd = 3;
        correlation{ii} = potTemp{2};
      else
        fprintf('Erroneous correlation, setting to mf \n')
        correlation{ii} = 'mf';
        strNum = 1;
        paramInd = 2;
      end
    else
      strNum = 1;
      paramInd = 2;
    end
    % take care of parameters
    potInputs(ii) = length(potTemp) - strNum;
    paramSingMat = potTemp{paramInd};
    for jj = paramInd+1:length(potTemp)
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
      if interactPotFlag
        runTemp{jj} = { names{jj}, correlation{jj}, paramSepCell{ii,jj} };
        strTemp = [  strTemp names{jj} correlation{jj} ...
          num2str( paramSepCell{ii,jj}, '_%.2f') '_' ];
      else
        runTemp{jj} = { names{jj}, paramSepCell{ii,jj} };
        strTemp = [  strTemp names{jj}  ...
          num2str( paramSepCell{ii,jj}, '_%.2f') '_' ];
      end
      % get rid of blanks just in case
      strTemp = strTemp( strTemp ~= ' ');
    end
    % fix str
    strTemp(end) = [];
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
