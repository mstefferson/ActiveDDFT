function plotOps( OPs, x, y, cScale, figStyle )
if nargin == 4
  figStyle = 'normal';
end
% set-up fonts
fontSize = 30;
% Find ticks
nx = length(x);
ny = length(y);
lx = x(end) - 2*x(1) + x(2);
ly = y(end) - 2*y(1) + y(2);
xMid = x( nx/2 + 1);
yMid = y( ny/2 + 1);
xTick = [xMid-lx/4 xMid xMid+lx/4];
yTick = [yMid-ly/4 yMid xMid+ly/4];
xLim  = [x(1) x(end)];
yLim  = [y(1) y(end)];
if nargin < 4
  cScale = 1;
end
% set-up subplot
numRow = 1;
numCol = 3;
% Set up a index vector so quiver is too crowded
divNumX = 8;
divNumY = 8;
deltaX  = ceil(nx / divNumX );
deltaY  = ceil(ny / divNumY);
% dir 1 = rows = x
subInd1 = 1:deltaX:(nx + 1 - deltaX);
% dir 2 = columns = y
subInd2 = 1:deltaY:(ny + 1 - deltaX);
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
% Concentration
axh1 = subplot(numRow,numCol,1); % Save the handle of the subplot
axh1.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh1);
h.TickLabelInterpreter = 'latex';
axh1.NextPlot = 'replaceChildren';
minC = min( OPs.C(:) );
maxC = max( OPs.C(:) );
if minC >= maxC -0.0001
  maxC = 1.1 .* minC ;
  minC = 0.9 .* minC;
end
axh1.CLim = cScale * [minC maxC];
axh1.YTick = yTick;
axh1.YLim = yLim;
axh1.XTick = xTick;
axh1.XLim = xLim;
wantedTickLabel =  num2cell( yTick ) ;
axh1.YTickLabel =  wantedTickLabel;
shading(axh1,'interp');
xlabel(axh1,'$$ x $$'); ylabel(axh1,'$$ y $$') 
axh1.FontSize = fontSize;
axis(axh1, 'square')
cTitle = '$$ C $$';
% Polar order
axh2 = subplot(numRow,numCol,2); % Save the handle of the subplot
axh2.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh2);
h.TickLabelInterpreter = 'latex';
axh2.NextPlot = 'replaceChildren';
axh2.CLim = [0 1];
axh2.XLim = xLim; %row and columns are flipped
axh2.YLim = yLim; %row and columns are flipped
axh2.XTick = xTick;
axh2.YTick = yTick;
axh2.YTickLabel =  wantedTickLabel;
shading(axh2,'interp');
xlabel(axh2,'$$ x $$'); ylabel(axh2,'$$ y $$') 
axh2.FontSize = fontSize;
axis(axh2, 'square')
pTitle = '$$ P $$';
% Nematic order
axh3 = subplot(numRow,numCol,3); % Save the handle of the subplot
axh3.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh3);
h.TickLabelInterpreter = 'latex';
axh3.NextPlot = 'replaceChildren';
axh3.CLim = [0 1];
axh3.XLim = xLim; %row and columns are flipped
axh3.YLim = yLim; %row and columns are flipped
axh3.XTick = xTick;
axh3.YTick = yTick;
axh3.YTickLabel = wantedTickLabel;
shading(axh3,'interp');
xlabel(axh3,'$$ x $$'); ylabel(axh3,'$$ y $$') 
axh3.FontSize = fontSize;
axis(axh3, 'square')
nTitle = '$$ N $$';
% Scale order parameters by it's max value to for it changes.
cTemp =  cScale * OPs.C;
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
