function plotBandSlice( cWant, fWant, bandTable )
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
wantedColors = getPlotLineColors( fd2Plot, 'log', 'viridis' );
% set-up figure
fig = figure();
fig.WindowStyle = 'normal';
fig.Position = [104 368 1015 263];
ax1 = subplot(1,3,1);
setAxBandSlice( ax1, xLabel, yLabels{1}, myTitle{1} );
ax2 = subplot(1,3,2);
setAxBandSlice( ax2, xLabel, yLabels{2}, myTitle{2} );
ax3 = subplot(1,3,3);
setAxBandSlice( ax3, xLabel, yLabels{3}, myTitle{3} );
%
for ii = 1:length(fd2Plot)
  %%
  % get max for shifting
  [~, maxInd] = max(cSlices2Plot{ii});
  circAmount = mod( n/2 - maxInd, n );
  % grab data
  c2plot = shiftSlice( cSlices2Plot{ii} ./ cWant, circAmount );
  p2plot = shiftSlice( pSlices2Plot{ii}, circAmount );
  n2plot = shiftSlice( nSlices2Plot{ii}, circAmount );
  keyboard
  % plot it
  plot(ax1, x, c2plot, 'Color', wantedColors(ii,:) )
  plot(ax2, x, p2plot, 'Color', wantedColors(ii,:) )
  plot(ax3, x, n2plot, 'Color', wantedColors(ii,:) )
end
legendCell = cellstr(num2str(fd2Plot, '$$Pe$$=%-d'));
hl = legend( legendCell, 'location', 'best' );
hl.Interpreter = 'latex';
hl.Position = [0.9013 0.3592 0.0872 0.3494];
%% functions
  function setAxBandSlice( ax, xLabel, yLabel, myTitle )
    fontSize = 16;
    ax.FontSize = fontSize;
    axis square
    xlabel(ax, xLabel)
    ylabel(yLabel)
    ax.TickLabelInterpreter = 'latex';
    axis(ax, 'square')
    hold(ax,'on')
    title(myTitle, 'Units', 'normalized', ...
      'Position', [0 1 0], 'HorizontalAlignment', 'left')
    box on
  end

  function data2plot = shiftSlice( data2plot, circAmount )
    % make sure everything is row vector
    if iscolumn(data2plot); data2plot = data2plot'; end
    % shift it to the center
    data2plot = circshift( data2plot, circAmount, 2 );
  end
end
