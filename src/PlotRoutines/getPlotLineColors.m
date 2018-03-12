% get colors based on nu/lplc values to indicate value
function wantedColors = getPlotLineColors( plotVec, scaleType, colorScheme )
if nargin == 2
  colorScheme = 'viridis';
end
% set up colors
fac = 255;
if strcmp( colorScheme, 'viridis' )
  colorArray = viridis( fac+1 );
elseif strcmp( colorScheme, 'pmkmp' )
  colorArray = pmkmp( fac+1 );
else
  fprintf('Do no recognize color scheme. Setting to viridis\n')
  colorArray = viridis( fac+1 );
end
% get size
lengthPv = length( plotVec );
% handle zeros separately
if strcmp( scaleType, 'log' )
  indsOk = plotVec > 0;
  plotVec = plotVec( indsOk);
  logVec = log10( plotVec );
  logVec = logVec - min( logVec );
  getInds = round( fac / max( logVec ) ...
    * ( logVec ) ) + 1;
  getInds( isinf(getInds) ) = fac+1;
elseif strcmp( scaleType, 'linear' )
  indsOk = ~isinf( plotVec );
  plotVec = plotVec( indsOk);
  getInds = round( fac / max( plotVec ) .* plotVec )+1;
else 
  error('Do not recognize the scale')
end
wantedColors = zeros( lengthPv, 3 );
wantedColors(indsOk,:) = colorArray( getInds, :);