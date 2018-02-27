function plotOps( OPs, x, y, cScale, figStyle )
if nargin == 4
  figStyle = 'normal';
end
% set-up grid things like ticks
[xTick, yTick, xLim, yLim, subInd1, subInd2] = buildTicks( x, y );
% Set up figure, make it a square 0.8 of
% smallest screen dimension
screenSize = get(0,'screensize');
screenWidth = screenSize(3); ScreenHeight = screenSize(4);
figWidth    = floor( screenWidth  );
figHeight   =  floor( ScreenHeight * .40);
figPos      = [ floor( 0.5 * ( screenWidth - figWidth ) ) ...
  floor( 0.5 * (ScreenHeight - figHeight ) ) ...
  figWidth figHeight];
%Build a square box set by smallest dimension of screen
Fig = figure();
Fig.WindowStyle = figStyle;
Fig.Position = figPos;
% Scale order parameters by it's max value to for it changes.
% set-up subplot
numRow = 1;
numCol = 3;
cTitle = 'C';
pTitle = 'P';
nTitle = 'N';
myTitle = {cTitle,pTitle,nTitle};
% Concentration
cTemp =  cScale * OPs.C;
axh1 = subplot(numRow,numCol,1); % Save the handle of the subplot
fixAxis( axh1, cTemp, xTick, yTick, xLim, yLim, myTitle{1} )
% Polar order
axh2 = subplot(numRow,numCol,2); % Save the handle of the subplot
fixAxis( axh2, [0 1], xTick, yTick, xLim, yLim, myTitle{2} )
% Nematic order
axh3 = subplot(numRow,numCol,3); % Save the handle of the subplot
fixAxis( axh3, [0 1], xTick, yTick, xLim, yLim, myTitle{3} )
% Scale order parameters by it's max value to for it changes.
polarTempX = OPs.POPx(subInd1,subInd2,:);
polarTempY = OPs.POPy(subInd1,subInd2,:);
nemTempX = OPs.NOPx(subInd1,subInd2,:);
nemTempY = OPs.NOPy(subInd1,subInd2,:);
nemTemp = OPs.NOP(subInd1,subInd2,:);
% scale nematic eigenvectors by their value
nemTempX = nemTempX .* nemTemp;
nemTempY = nemTempY .* nemTemp;
% These matrices will need to be transposed to correct for x and y
% loop over frames
ii = 1;
%% Concentration
subplot(axh1);
cla(axh1);
pcolor( axh1, x, y, cTemp(:,:,ii)' );
shading interp
TitlStr = cTitle;
title(axh1,TitlStr);
pause(0.001);
drawnow;
%% Polar order
subplot(axh2);
cla(axh2);
pcolor(axh2, x, y, OPs.POP(:,:,ii)' );
shading interp
TitlStr = pTitle;
hold on
quiver(axh2, x(subInd1), y(subInd2),...
  polarTempX(:,:,ii)',...
  polarTempY(:,:,ii)', 0,'color',[1,1,1] );
title(axh2,TitlStr)
%% Nematic order
cla(axh3)
subplot(axh3);
pcolor(axh3, x, y, OPs.NOP(:,:,ii)');
shading interp
TitlStr = nTitle;
hold on
quiver(axh3, x(subInd1),y(subInd2),...
  nemTempX(:,:,ii)', nemTempY(:,:,ii)',0,...
  'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
title(axh3,TitlStr)
end
% set-up grid things like ticks
function [xTick, yTick, xLim, yLim, subInd1, subInd2] = buildTicks( x, y )
nx = length(x);
ny = length(y);
lx = x(end) - 2*x(1) + x(2);
ly = y(end) - 2*y(1) + y(2);
% Find ticks
xMid = x( nx/2 + 1);
yMid = y( ny/2 + 1);
xTick = [xMid-lx/4 xMid xMid+lx/4];
yTick = [yMid-ly/4 yMid xMid+ly/4];
xLim  = [x(1) x(end)];
yLim  = [y(1) y(end)];
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
% fix axis
function fixAxis( ax, data, xTick, yTick, xLim, yLim, myTitle )
fontSize = 30;
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
ax.CLim = [minVal maxVal];
ax.YTick = yTick;
ax.YLim = yLim;
ax.XTick = xTick;
ax.XLim = xLim;
wantedTickLabel =  num2cell( yTick ) ;
ax.YTickLabel =  wantedTickLabel;
shading(ax,'interp');
xlabel(ax,'$$ x $$'); ylabel(ax,'$$ y $$')
title(myTitle)
axis(ax, 'square')
end

