function noiseObj = noiseRunManager( noiseCell, n3 )
% no rotational noise if 2D
if n3 == 1
  noiseCell{2} = 0;
end
% build noise obj
if [noiseCell{1} noiseCell{2} ] == 0
  noiseInd = 1;
  noiseStr = {''};
  noiseRuns = { [0 0] };
else
  noiseMat = combvec( noiseCell{1}, noiseCell{2} )';
  numNoise = size( noiseMat, 1 );
  noiseInd = 1:numNoise;
  noiseStr = cell(1,numNoise);
  noiseRuns = cell(1,numNoise);
  for ii = 1:numNoise
    noiseStr{ii} = ['_noisep' ...
      num2str( noiseMat(ii,1), '%.1f') 'r' num2str(noiseMat(ii,2), '%.1f') ];
    noiseRuns{ii} = noiseMat(ii,:);
  end
end
noiseObj.param = noiseRuns;
noiseObj.str = noiseStr;
noiseObj.inds = noiseInd;
