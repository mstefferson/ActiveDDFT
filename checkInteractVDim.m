% check potential
function pot = checkInteractVDim( pot, n3 )
for ii = 1:length(pot)
  currPot = pot{ii};
  if strcmp( currPot{1}, 'pa2d' ) && n3 == 1
    fprintf('Deleting potential %s because no n3\n',...
      currPot{1});
    pot{ii} = [];
  end
end
pot = pot( ~cellfun('isempty', pot) );
