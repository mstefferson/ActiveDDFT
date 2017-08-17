% check potential
function pot = checkExtVDim( pot, nVec )
for ii = 1:length(pot)
  currPot = pot{ii};
  if strcmp(currPot{1},'nemV') 
    if nVec(3) == 1
      fprintf('Deleting potential %s. No 3 dimension \n',...
      currPot{1} );
      pot{ii} = [];
    end
  elseif strcmp(currPot{1},'polV') 
    if nVec(3) == 1
      fprintf('Deleting potential %s. No 3 dimension\n',...
      currPot{1} );
      pot{ii} = [];
    end
  elseif nVec( currPot{2}(1) ) == 1
    fprintf('Deleting potential %s. No dimension %d\n',...
      currPot{1}, currPot{2}(1) );
    pot{ii} = [];
  end
end
pot = pot( ~cellfun('isempty', pot) );
