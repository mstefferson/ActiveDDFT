% Makes movie of OP vs time
function OPMovieMakerTgtherDirAvi(MovStr,x,y,...
  OP,TimeRec, cScale, rhoSlice)
% set font-size
fontSize = 30;
% set colormap
randFig = randi(10000);
figure(randFig);
colormap( viridis );
close(randFig);
% set up dots
cMat = copper(3);
% frames and things
nFrames = length(TimeRec);
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
% Set up figure, make it a square 0.8 of
% smallest screen dimension
screenSize = get(0,'screensize');
screenWidth = screenSize(3); ScreenHeight = screenSize(4);
figWidth    = floor( screenWidth );
figHeight   =  floor( ScreenHeight * .40);
figPos      = [ floor( 0.5 * ( screenWidth - figWidth ) ) ...
  floor( 0.5 * (ScreenHeight - figHeight ) ) ...
  figWidth figHeight];
%Build a square box set by smallest dimension of screen
Fig = figure();
Fig.WindowStyle = 'normal';
Fig.Position = figPos;
%%
%Initialize the movie structure
Mov = VideoWriter(MovStr);
Mov.FrameRate = 4;
open(Mov);
% set-up subplot
numRow = 1;
numCol = 3;
% Concentration
set(gcf,'renderer','zbuffer')
axh1 = subplot(numRow,numCol,1); % Save the handle of the subplot
axh1.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh1);
h.TickLabelInterpreter = 'latex';
axh1.NextPlot = 'replaceChildren';
minC = min( OP.C_rec(:) );
maxC = max( OP.C_rec(:) );
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
axis square
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
%axh3.YDir = 'rev';
axh3.XTick = xTick;
axh3.YTick = yTick;
axh3.YTickLabel = wantedTickLabel;
shading(axh3,'interp');
xlabel(axh3,'$$ x $$'); ylabel(axh3,'$$ y $$') 
axh3.FontSize = fontSize;
axis square
nTitle = '$$ N $$';
% inset
if rhoSlice.plotInset
  axhinsetMain = axh1;
  % set-up inset
  posVec = axhinsetMain.Position;
  x0 = posVec(1);
  y0 = posVec(2);
  w0 = posVec(3);
  h0 = posVec(4);
  axLabColor = [0.85 0 0];
  % Place second set of axes on same plot
  x1 = x0 + w0 / (1.80);
  y1 = y0 + h0 / (2.1);
  w1 = w0 / (2.5);
  h1 = h0 / (2.0);
  axhinset = axes('Position', [x1 y1 w1 h1]);
  axis square
  posVec = axhinset.Position;
  axhinset.Position = posVec ;
  axhinset.XColor = axLabColor;
  axhinset.YColor = axLabColor;
  axhinset.XLabel.Color = axLabColor;
  axhinset.YLabel.Color = axLabColor;
  axhinset.YLabel.String = '$$ f ( \phi ) $$';
  axhinset.XLabel.String = '$$ \phi $$';
  yMinInset = ...
    min( [ rhoSlice.slice1(:); rhoSlice.slice2(:); rhoSlice.slice3(:) ] );
  yMaxInset = ...
    max( [ rhoSlice.slice1(:); rhoSlice.slice2(:); rhoSlice.slice3(:) ] );
  axhinset.YLim = [yMinInset yMaxInset];
  axhinset.NextPlot = 'replaceChildren';
end
% Scale order parameters by it's max value to for it changes.
cTemp =  cScale * OP.C_rec;
maxP = max( OP.POP_rec(:) );
polarTempX = OP.POPx_rec(subInd1,subInd2,:) ./ maxP;
polarTempY = OP.POPy_rec(subInd1,subInd2,:) ./ maxP;
nemTempX = OP.NOPx_rec(subInd1,subInd2,:);
nemTempY = OP.NOPy_rec(subInd1,subInd2,:);
nemTemp = OP.NOP_rec(subInd1,subInd2,:);
% scale nematic eigenvectors by their value
nemTempX = nemTempX .* nemTemp;
nemTempY = nemTempY .* nemTemp;
% These matrices will need to be transposed to correct for x and y
% loop over frames
try
  vec2loop = 1:nFrames;
  for ii = vec2loop
    % get t in title
    tTitle = [' ($$ t = $$ ' num2str( TimeRec(ii), '%.2f' ) ')' ];
    %% Concentration
    subplot(axh1);
    cla(axh1);
    pcolor( axh1, x, y, cTemp(:,:,ii)' );
    shading interp
    TitlStr = [cTitle  tTitle]; 
    title(axh1,TitlStr);
    pause(0.001);
    drawnow;
    %% Polar order
    subplot(axh2);
    cla(axh2);
    pcolor(axh2, x, y, OP.POP_rec(:,:,ii)' );
    shading interp
    TitlStr = [pTitle  tTitle];  
    hold on
    quiver(axh2, x(subInd1), y(subInd2),...
      polarTempX(:,:,ii)',...
      polarTempY(:,:,ii)', 0,'color',[1,1,1] );
    title(axh2,TitlStr)
    %% inset
    if rhoSlice.plotInset
      % x's
      hold(axhinsetMain,'on')
      p = plot( axhinsetMain, x(rhoSlice.dotx1), y(rhoSlice.doty1),'x','MarkerSize',20 );
      p.Color = cMat(1,:);
      p = plot( axhinsetMain, x(rhoSlice.dotx2), y(rhoSlice.doty2),'x','MarkerSize',20 );
      p.Color = cMat(2,:);
      p = plot( axhinsetMain, x(rhoSlice.dotx3), y(rhoSlice.doty3),'x','MarkerSize',20 );
      p.Color = cMat(3,:);
      % distribution inset
      cla(axhinset);
      p = plot( axhinset, rhoSlice.phi, rhoSlice.slice1(:,ii) );
      hold( axhinset , 'on' )
      p.Color = cMat(1,:);
      p = plot( axhinset, rhoSlice.phi, rhoSlice.slice2(:,ii) );
      p.Color = cMat(2,:);
      p = plot( axhinset, rhoSlice.phi, rhoSlice.slice3(:,ii) );
      p.Color = cMat(3,:);
      hold( axhinset, 'off')
      axhinset.XLim = [ rhoSlice.phi(1) rhoSlice.phi(end) ];
      axhinset.Position = posVec ;
      axhinset.XColor = axLabColor;
      axhinset.YColor = axLabColor;
      axhinset.XLabel.Color = axLabColor;
      axhinset.YLabel.Color = axLabColor;
      axhinset.YLabel.String = '$$ f ( \phi ) $$';
      axhinset.XLabel.String = '$$ \phi $$';
      axis( axhinset, 'square' )
      axhinset.YLim = [yMinInset yMaxInset];
    end
    pause(0.01);
    drawnow;
    %% Nematic order
    cla(axh3)
    subplot(axh3);
    pcolor(axh3, x, y, OP.NOP_rec(:,:,ii)');
    shading interp
    TitlStr = [nTitle  tTitle]; 
    hold on
    quiver(axh3, x(subInd1),y(subInd2),...
      nemTempX(:,:,ii)', nemTempY(:,:,ii)',0,...
      'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
    title(axh3,TitlStr)
    pause(0.001);
    drawnow;
    % Fig.Position can drift. Has to do with capturing a figure
    % before graphics render. Fix with drawnow and a pause.
    Fr = getframe(Fig);
    writeVideo(Mov,Fr);
  end% End frame loop
catch err
  fprintf('%s', err.getReport('extended')) ;
end % try catch
close(Fig)
