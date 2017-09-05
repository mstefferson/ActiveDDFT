function opGetSnapshot( path2dir, frame )
% set font-size
fontSize = 30;
% set-up files
opName = dir( [path2dir '/op_*'] );
runName = dir( [path2dir '/run_*'] );
paramName = dir( [path2dir '/params_*'] );
oPFile = matfile([opName.folder '/' opName.name]);
load( [paramName.folder '/' paramName.name] )
runFile = matfile([runName.folder '/' runName.name]);
gridObj = runFile.gridObj;
cScale = particleObj.b;
% store files in ram
oP.C_rec = oPFile.C_rec;
oP.POP_rec = oPFile.POP_rec;
oP.POPx_rec = oPFile.POPx_rec;
oP.POPy_rec = oPFile.POPy_rec;
oP.NOP_rec = oPFile.NOP_rec;
oP.NOPx_rec = oPFile.NOPx_rec;
oP.NOPy_rec = oPFile.NOPy_rec;
if isfield(oP, 'sliceRho')
  rhoSlice = oP.sliceRho;
else
  rhoSlice.plotInset = 0;
end
% fix frame if too large
maxFrame = length( oPFile.OpTimeRecVec );
frame = min( frame, maxFrame );
%
x = gridObj.x1;
y = gridObj.x2;
% set-up fig
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
DivNumX = 8;
DivNumY = 8;
DeltaX  = ceil(nx / DivNumX );
DeltaY  = ceil(ny / DivNumY);
% dir 1 = rows = x
SubInd1 = 1:DeltaX:(nx + 1 - DeltaX);
% dir 2 = columns = y
SubInd2 = 1:DeltaY:(ny + 1 - DeltaX);
% Set up figure, make it a square 0.8 of
% smallest screen dimension
ScreenSize = get(0,'screensize');
ScreenWidth = ScreenSize(3); ScreenHeight = ScreenSize(4);
FigWidth    = floor( ScreenWidth  );
FigHeight   =  floor( ScreenHeight * .40);
FigPos      = [ floor( 0.5 * ( ScreenWidth - FigWidth ) ) ...
  floor( 0.5 * (ScreenHeight - FigHeight ) ) ...
  FigWidth FigHeight];
%Build a square box set by smallest dimension of screen
Fig = figure();
Fig.WindowStyle = 'normal';
Fig.Position = FigPos;
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
minC = min( oP.C_rec(:) );
maxC = max( oP.C_rec(:) );
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
%axh3.YDir = 'rev';
axh3.XTick = xTick;
axh3.YTick = yTick;
axh3.YTickLabel = wantedTickLabel;
shading(axh3,'interp');
xlabel(axh3,'$$ x $$'); ylabel(axh3,'$$ y $$') 
axh3.FontSize = fontSize;
axis(axh3, 'square')
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
cTemp =  cScale * oP.C_rec;
polarTempX = oP.POPx_rec(SubInd1,SubInd2,:);
polarTempY = oP.POPy_rec(SubInd1,SubInd2,:);
maxPolar = max( max( max( oP.POP_rec(SubInd1,SubInd2,:) ) ) );
polarTempX = polarTempX ./ maxPolar;
polarTempY = polarTempY ./ maxPolar;
nemTempX = oP.NOPx_rec(SubInd1,SubInd2,:);
nemTempY = oP.NOPy_rec(SubInd1,SubInd2,:);
maxNem = max( max( max( oP.NOP_rec(SubInd1,SubInd2,:) ) ) );
nemTempX = nemTempX .* oP.NOP_rec(SubInd1,SubInd2,:) ./ maxNem;
nemTempY = nemTempY .* oP.NOP_rec(SubInd1,SubInd2,:) ./ maxNem;
% These matrices will need to be transposed to correct for x and y
% loop over frames
ii = frame;
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
pcolor(axh2, x, y, oP.POP_rec(:,:,ii)' );
shading interp
TitlStr = pTitle;
hold on
quiver(axh2, x(SubInd1), y(SubInd2),...
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
%% Nematic order
cla(axh3)
subplot(axh3);
pcolor(axh3, x, y, oP.NOP_rec(:,:,ii)');
shading interp
TitlStr = nTitle;
hold on
quiver(axh3, x(SubInd1),y(SubInd2),...
  nemTempX(:,:,ii)', nemTempY(:,:,ii)',0,...
  'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
title(axh3,TitlStr)
