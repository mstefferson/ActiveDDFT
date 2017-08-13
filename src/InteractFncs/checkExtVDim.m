% check potential
function pot = checkExtVDim( pot, nVec )
for ii = 1:length(pot)
  currPot = pot{ii};
  if nVec( currPot{2}(1) ) == 1
    fprintf('Deleting potential %s wrong dimension %d\n',...
      currPot{1}, currPot{2}(1) );
    pot{ii} = [];
  end
end
pot = pot( ~cellfun('isempty', pot) );
