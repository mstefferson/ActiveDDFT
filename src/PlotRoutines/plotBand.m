function plotBand(fileName )
% titles
myTitle = {'A', 'B', 'C'};
% set some figur parameters
numRow = 1;
numCol = 3;
% load params and opObj from current WD
paramsMat = dir( [ fileName '/params_*' ] );
opMat = dir( [ fileName '/op_*' ] );
paramLoad = load( [ fileName '/' paramsMat.name ] );
fileObj = matfile( [ fileName   '/' opMat.name ] );
% scale c
b = 1 / pi;
fprintf('Scaling C by b = %f\n', b );
C = fileObj.C_rec(:,:,end) * b;
P = fileObj.POP_rec(:,:,end);
N = fileObj.NOP_rec(:,:,end);
% set-up figure
fig = figure();
fig.WindowStyle = 'normal';
fig.Position = [2 361 1233 298];
% set-up grid things like ticks
[xTick, yTick, xLim, yLim, subInd1, subInd2, x1, x2] = ....
  buildTicks( paramLoad.systemObj );
% Concentration
axh1 = subplot(numRow,numCol,1); % Save the handle of the subplot
fixAxis( axh1, C, xTick, yTick, xLim, yLim, myTitle{1} )
% Polar order
axh2 = subplot(numRow,numCol,2); % Save the handle of the subplot
fixAxis( axh2, P, xTick, yTick, xLim, yLim, myTitle{2} )
% Nematic order
axh3 = subplot(numRow,numCol,3); % Save the handle of the subplot
fixAxis( axh3, N, xTick, yTick, xLim, yLim, myTitle{3} )
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

  function fixAxis( ax, data, xTick, yTick, xLim, yLim, myTitle )
    fontSize = 16;
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
    title(myTitle, 'Units', 'normalized', ...
      'Position', [0 1 0], 'HorizontalAlignment', 'left')
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
end
