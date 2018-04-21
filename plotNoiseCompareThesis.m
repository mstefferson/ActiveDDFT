%% Plot of external potential test for thesis
dirname = '/Users/mike/Projects/ActiveDDFT/analyzedfiles/Noise4ThesisPlots';
my_dirs = dir([dirname '/Hr*']);
saveName = 'noise_compare';
saveMe = 1;
%%
% titles
fig = figure();
fig.WindowStyle = 'normal';
fig.Position = [472 -206 1036 646];
myTitle = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'};
% set some figur parameters
numRow = 3;
numCol = 3;
for ii = 1:numRow
  plotRowId = numRow*(ii-1);
  % load params and opObj from current WD
  dir_temp = [my_dirs(ii).folder '/' my_dirs(ii).name];
  fileNameTemp = [dir_temp '/'];
  fileName = dir(fileNameTemp);
  paramsMat = dir( [ fileNameTemp '/params_*' ] );
  opMat = dir( [ fileNameTemp '/op_*' ] );
  paramLoad = load( [ fileNameTemp '/' paramsMat.name ] );
  fileObj = matfile( [ fileNameTemp   '/' opMat.name ] );
  % scale c
  b = 1 / pi;
  fprintf('Scaling C by b = %f\n', b );
  C = fileObj.C_rec(:,:,end) * b;
  P = fileObj.POP_rec(:,:,end);
  N = fileObj.NOP_rec(:,:,end);
  % set-up figure
  % set-up grid things like ticks
  [xTick, yTick, xLim, yLim, subInd1, subInd2, x1, x2] = ....
    buildTicks( paramLoad.systemObj );
  % Concentration
  axh1 = subplot(numRow,numCol,plotRowId+1); % Save the handle of the subplot
  % super hack axis
  %   cLim = [0.6 3.8];
  cLim = [min(C(:)) max(C(:))];
  fixAxis( axh1, C, xTick, yTick, xLim, yLim, cLim, myTitle{plotRowId+1} )
  % Polar order
  axh2 = subplot(numRow,numCol,plotRowId+2); % Save the handle of the subplot
  % super hack axis
  %   cLim = [0 0.1]
  if ii ~= 3
    cLim = [0 2e-3];
  else
    cLim = [min(P(:)) max(P(:))];
  end
  fixAxis( axh2, P, xTick, yTick, xLim, yLim, cLim, myTitle{plotRowId+2} )
  % Nematic order
  if ii == 1
    cLim = [0 2e-3];
  else
    cLim = [min(N(:)) max(N(:))];
  end
  axh3 = subplot(numRow,numCol,plotRowId+3); % Save the handle of the subplot
  fixAxis( axh3, N, xTick, yTick, xLim, yLim, cLim, myTitle{plotRowId+3} )
  % shift things
  [circAmount, shiftDim] = findShift( C, paramLoad.systemObj );
  C = shiftData( C, circAmount, shiftDim );
  P = shiftData( P, circAmount, shiftDim );
  N = shiftData( N, circAmount, shiftDim );
  pX = shiftData( fileObj.POPx_rec(:,:,end), circAmount, shiftDim );
  pY = shiftData( fileObj.POPy_rec(:,:,end), circAmount, shiftDim );
  nX = shiftData( fileObj.NOPx_rec(:,:,end), circAmount, shiftDim );
  nY = shiftData( fileObj.NOPy_rec(:,:,end), circAmount, shiftDim );
  % set up arrows
  maxP = max( P(:) );
  polarTempX = pX(subInd1,subInd2,end) ./ maxP;
  polarTempY = pY(subInd1,subInd2,end) ./ maxP;
  nemTempX = nX(subInd1,subInd2,end);
  nemTempY = nY(subInd1,subInd2,end);
  nemTemp = N(subInd1,subInd2,end);
  % scale nematic eigenvectors by their value
  nemTempX = nemTempX .* nemTemp;
  nemTempY = nemTempY .* nemTemp;
  % These matrices will need to be transposed to correct for x and y
  %% Concentration
  pcolor( axh1, x1, x2, C' );
  shading( axh1, 'interp')
  %% Polar order
  pcolor(axh2, x1, x2, P' );
  shading( axh2, 'interp')
  pause(1)
  drawnow
  %% Nematic order
  pcolor(axh3, x1, x2, N.');
  shading( axh3, 'interp')
  pause(0.1)
  drawnow
  quiver(axh2, x1(subInd1), x2(subInd2),...
    polarTempX', polarTempY', 0,'color',[1,1,1] );
  quiver(axh3, x1(subInd1),x2(subInd2),...
    nemTempX', nemTempY',0,...
    'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
end
if saveMe
  saveas(gcf, saveName,'png')
end
%% functions
function [xTick, yTick, xLim, yLim, subInd1, subInd2, x1, x2] = ...
  buildTicks( systemObj )
% build grid
nx = systemObj.n1;
ny = systemObj.n2;
x1 = systemObj.l1 / nx .* ( -nx/2:(nx/2-1) );
x2 = systemObj.l2 / ny .* ( -ny/2:(ny/2-1) );
lx = systemObj.l1;
ly = systemObj.l2;
% Find ticks
xMid = x1( nx/2 + 1);
yMid = x2( ny/2 + 1);
xTick = [xMid-lx/4 xMid xMid+lx/4];
yTick = [yMid-ly/4 yMid xMid+ly/4];
xLim  = [x1(1) x1(end)];
yLim  = [x2(1) x2(end)];
% Set up a index vector so quiver is too crowded
divNumX = 8;
divNumY = 8;
deltaX  = ceil(nx / divNumX );
deltaY  = ceil(ny / divNumY);
% dir 1 = rows = x
subInd1 = 1:deltaX:(nx + 1 - deltaX);
% dir 2 = columns = y
subInd2 = 1:deltaY:(ny + 1 - deltaX);
end

function fixAxis( ax, data, xTick, yTick, xLim, yLim, cLim, myTitle )
fontSize = 12;
ax.NextPlot = 'replaceChildren';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = fontSize;
hold(ax, 'on')
h = colorbar('peer',ax);
h.TickLabelInterpreter = 'latex';
minVal = min( data(:) );
maxVal = max( data(:) );
if minVal >= maxVal -0.0001
  maxVal = 1.1 .* minVal ;
  minVal = 0.9 .* minVal;
end
%ax.CLim = [minVal maxVal];
ax.CLim = cLim;
ax.YTick = yTick;
ax.YLim = yLim;
ax.XTick = xTick;
ax.XLim = xLim;
wantedTickLabel =  num2cell( yTick ) ;
ax.YTickLabel =  wantedTickLabel;
shading(ax,'interp');
xlabel(ax,'$$ x $$'); ylabel(ax,'$$ y $$')
ht = title(ax, myTitle);
% hack in title position
ht.Position(1) = ht.Position(1) - 9*ht.Position(3);
ht.Position(2) = ht.Position(2) - 0.3;
ht.FontSize = 14;
% for some reason, shifting the title to the left fucks
% up the title. it disappears
% ht = title(ax, myTitle, 'Units', 'normalized', ...
%   'Position', [0 1 0], 'HorizontalAlignment', 'left');
axis(ax, 'square')
end
function [circAmount, shiftDim] = findShift( dataMap, systemObj )
% shift things
[~, rows, cols, ~, ~, ~] = ...
  sliceInhomoVar( systemObj, dataMap );
if length(rows) > length(cols)
  shiftDim =  1;
  [~,maxInd] = max( dataMap(:,1) );
  circAmount = mod( systemObj.n1/2 - maxInd, systemObj.n1 );
else
  shiftDim = 2;
  [~,maxInd] = max( dataMap(1,:) );
  circAmount = mod( systemObj.n2/2 - maxInd, systemObj.n2 );
end
end
function data2plot = shiftData( data2plot, circAmount, shiftDim )
% shift it to the center
data2plot = circshift( data2plot, circAmount, shiftDim );
end

