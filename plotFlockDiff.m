dirname = '/Users/mike/Projects/ActiveDDFT/analyzedfiles/FlockDiff';
my_dirs = dir([dirname '/Hr*']);
saveName = 'flock_diff_break_apart';
pLim = [0 0.07];
% saveName = 'flock_diff_backscatter';
% pLim = [0 0.175];
saveMe = 1;
% myTitle = {'A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'D1', 'D2',...
%   'E1', 'E2', 'F1', 'F2', 'G1', 'G2', 'H1', 'H2'};
myTitle = {'C0',  'C1.1', 'C1.2', 'P0','P1.1', 'P1.2',...
  'C0',  'C2.1', 'C2.2', 'P0', 'P2.1',  'P2.2'};
% myTitle= cell(100,1);
% time inds
timeId = [1 4 7];
% titles
% fig.WindowStyle = 'normal';
% set some figur parameters
numRow = 2*numT;
numCol = length(timeId);
numT = length(my_dirs);
fig = figure();
fig.WindowStyle = 'normal';
fig.Position = [426 1 754 704];
for ii = 1:numT
  % load params and opObj from current WD
  dir_temp = [my_dirs(ii).folder '/' my_dirs(ii).name];
  fileNameTemp = [dir_temp '/'];
  disp(my_dirs(ii).name)
  paramsMat = dir( [ fileNameTemp '/params_*' ] );
  opMat = dir( [ fileNameTemp '/op_*' ] );
  paramLoad = load( [ fileNameTemp '/' paramsMat.name ] );
  fileObj = matfile( [ fileNameTemp   '/' opMat.name ] );
  % scale c
  b = 1 / pi;
  fprintf('Scaling C by b = %f\n', b );
  C = fileObj.C_rec(:,:,:) * b;
  P = fileObj.POP_rec(:,:,:) * b;
  % set-up figure
  [xTick, yTick, xLim, yLim, subInd1, subInd2, x1, x2] = ....
    buildTicks( paramLoad.systemObj );
  nPlot = length(timeId);
  myColors = viridis(nPlot);
  maxP = max( P(:) );
  for nn = 1:nPlot
    plotid = timeId(nn);
    cTemp = C(:,:,plotid);
    pTemp = P(:,:,plotid);
    pX = fileObj.POPx_rec(:,:,plotid);
    pY = fileObj.POPy_rec(:,:,plotid);
    pXYnorm = 3*maxP;
    polarTempX = pX(subInd1,subInd2) ./ pXYnorm;
    polarTempY = pY(subInd1,subInd2) ./ pXYnorm;
    % Concentration
    plotidtemp = 2*(ii-1)*numCol+nn;
    %     disp(plotidtemp)
    axh1 = subplot(numRow,numCol,plotidtemp); % Save the handle of the subplot
    cLim = [min(C(:)) max(C(:))];
    fixAxis( axh1, C, xTick, yTick, xLim, yLim, cLim, ...
      myTitle{plotidtemp} )
    pcolor( axh1, x1, x2, cTemp' );
    shading( axh1, 'interp')
    %% Polar order
    plotidtemp = 2*(ii-1)*numCol+numCol+nn;
    %     disp(plotidtemp)
    axh2 = subplot(numRow,numCol,plotidtemp); % Save the handle of the subplot
    fixAxis( axh2, C, xTick, yTick, xLim, yLim, pLim, ...
      myTitle{plotidtemp})
    pcolor(axh2, x1, x2, pTemp' );
    shading( axh2, 'interp')
    pause(1)
    drawnow
    quiver(axh2, x1(subInd1), x2(subInd2),...
      polarTempX', polarTempY', 0,'color',[1,1,1] );
  end
end
if saveMe
%   disp([saveName num2str(ii)])
  saveas(fig, saveName,'png')
end
%%%%%%%%% functions %%%%%%%%%%%%%%%
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
% minVal = min( data(:) );
% maxVal = max( data(:) );
% if minVal >= maxVal -0.0001
%   maxVal = 1.1 .* minVal ;
%   minVal = 0.9 .* minVal;
% end
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
ht.Position(1) = ht.Position(1) - 8*ht.Position(3);
ht.Position(2) = ht.Position(2) - 0.3;
ht.FontSize = 12;
% for some reason, shifting the title to the left fucks
% up the title. it disappears
% ht = title(ax, myTitle, 'Units', 'normalized', ...
%   'Position', [0 1 0], 'HorizontalAlignment', 'left');
axis(ax, 'square')
end

