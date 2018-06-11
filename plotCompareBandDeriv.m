function plotCompareBandDeriv( cWant, fWant, bandTable )
% titles
myTitle = {'A', 'B', 'C'};
yLabels = { 'Concentration $$ C^* $$', 'Polar order $$ P $$', ...
  'Nematic order $$ N $$'};
xLabel = 'y';
% band slices
% get data
fd2Plot = bandTable.fd( bandTable.c == cWant ...
  & any(bandTable.fd == fWant,2) );
cSlices2Plot = bandTable.cSlice( bandTable.c == cWant...
  & any(bandTable.fd == fWant,2));
pSlices2Plot = bandTable.pSlice( bandTable.c == cWant...
  & any(bandTable.fd == fWant,2));
nSlices2Plot = bandTable.nSlice( bandTable.c == cWant...
  & any(bandTable.fd == fWant,2));
% sort it
[fd2Plot, indSort] = sort( fd2Plot );
cSlices2Plot = cSlices2Plot(indSort);
pSlices2Plot = pSlices2Plot(indSort);
nSlices2Plot = nSlices2Plot(indSort);
n = length( cSlices2Plot{1} );
l = 10;
x = l / n * (-n/2:n/2-1);
% get colors
% wantedColors = getPlotLineColors( fd2Plot, 'log', 'viridis' );
% build k
k = pi/l * [-n/2:(n/2-1)];
%
for ii = 1:length(fd2Plot)
  %%
  % set-up figure
  fig = figure();
  fig.WindowStyle = 'docked';
  ax1 = subplot(1,3,1);
  hold on
%   setAxBandSlice( ax1, xLabel, yLabels{1}, myTitle{1} );
  ax2 = subplot(1,3,2);
  hold on
%   setAxBandSlice( ax2, xLabel, yLabels{2}, myTitle{2} );
  ax3 = subplot(1,3,3);
  hold on
%   setAxBandSlice( ax3, xLabel, yLabels{3}, myTitle{3} );
  % get max for shifting
  [~, maxInd] = max(cSlices2Plot{ii});
  circAmount = mod( n/2 - maxInd, n );
  % grab data
  c2plot = shiftSlice( cSlices2Plot{ii} ./ cWant, circAmount );
  p2plot = shiftSlice( pSlices2Plot{ii}, circAmount );
  n2plot = shiftSlice( nSlices2Plot{ii}, circAmount );
  ck = fftshift(fftn(c2plot));
  dc = real(ifftn(ifftshift(sqrt(-1) * k .* ck)));
  % flip it
  dc(n/2:end) = -dc(n/2:end); 
  nk = fftshift(fftn(n2plot));
  dn = real(ifftn(ifftshift(sqrt(-1) * k .* nk)));
  dn(n/2:end) = -dn(n/2:end);
  % normalize
  c2plot = c2plot ./max(c2plot);
  dc = dc./max(dc);
  p2plot = p2plot ./max(p2plot);
  n2plot = n2plot ./max(n2plot);
  dn = dn./max(dn);
  % plot it
  %   plotyy(ax1, x, c2plot, x, dc, 'Color', wantedColors(ii,:) )
  plotyy(ax1, x, c2plot, x, dc)
  plot(ax2, x, dc,'r' )
  plot(ax2, x, p2plot, 'b')
%    plot(ax2, x, p2plot, 'Color', wantedColors(ii,:) )
  plot(ax2, x, dn,'k' )
  legend(ax2,'dc','p','dn', 'location', 'best')
  %   plotyy(ax3, x, n2plot, x, dn, 'Color', wantedColors(ii,:) )
  plotyy(ax3, x, n2plot./max(n2plot), x, dn./max(dn))
  title(num2str(ii))
end
% legendCell = cellstr(num2str(fd2Plot, '$$Pe$$=%-d'));
% hl = legend( legendCell, 'location', 'best' );
% hl.Interpreter = 'latex';
% hl.Position = [0.9013 0.3592 0.0872 0.3494];
%% functions
  function data2plot = shiftSlice( data2plot, circAmount )
    % make sure everything is row vector
    if iscolumn(data2plot); data2plot = data2plot'; end
    % shift it to the center
    data2plot = circshift( data2plot, circAmount, 2 );
  end
end
