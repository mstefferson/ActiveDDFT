% Makes movie of OP vs time
function OPMovieMakerTgtherDirAviFix(MovStr,x,y,phi,...
  OP,DistRec,TimeRec, cScale)
% matlab calls rows y and I call rows x. Flip mine to malabs
xT = x;
x = y;
y = xT;
% set colormap
colormap( viridis );
% frames and things
nFrames = length(TimeRec);
nx = length(x);
ny = length(y);
l1 = x(end) - 2*x(1) + x(2);
l2 = y(end) - 2*y(1) + y(2);
% Find ticks
xTick = round( [-l1/4 0 l1/4] );
yTick = round( [-l2/4 0 l2/4] );
xLim  = [x(1) x(end)];
yLim  = [y(1) y(end)];
% xLim = [-l1/2 l1/2];
% yLim = [-l2/2 l2/2];
% Set up a index vector so quiver is too crowded
DivNumX = 8;
DivNumY = 8;
DeltaX  = ceil(nx / DivNumX );
DeltaY  = ceil(ny / DivNumY);
% dir 1 = rows = y
SubInd1 = 1:DeltaY:(ny + 1 - DeltaY);
% dir 2 = columns = x
SubInd2 = 1:DeltaX:(nx + 1 - DeltaX);
% Set up figure, make it a square 0.8 of
% smallest screen dimension
ScreenSize = get(0,'screensize');
ScreenWidth = ScreenSize(3); ScreenHeight = ScreenSize(4);
FigWidth    = floor( ScreenWidth * .6 );
FigHeight   =  floor( ScreenHeight * .8);
FigPos      = [ floor( 0.5 * ( ScreenWidth - FigWidth ) ) ...
  floor( 0.5 * (ScreenHeight - FigHeight ) ) ...
  FigWidth FigHeight];
%Build a square box set by smallest dimension of screen
Fig = figure();
Fig.WindowStyle = 'normal';
Fig.Position = FigPos;
%Initialize the movie structure
Mov = VideoWriter(MovStr);
Mov.FrameRate = 4;
open(Mov);
% Concentration
set(gcf,'renderer','zbuffer')
axh1 = subplot(2,2,1); % Save the handle of the subplot
axh1.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh1);
h.TickLabelInterpreter = 'latex';
axh1.NextPlot = 'replaceChildren';
minC = min(min(min(OP.C_rec )));
maxC = max(max(max(OP.C_rec)));
if minC >= maxC
  minC = 0;
  maxC = 0.1;
end
axh1.CLim = cScale * [minC maxC];
% axh1.XLim = xLim; %row and columns are flipped
% axh1.YLim = yLim; %row and columns are flipped
%axh1.YDir = 'rev';
axh1.YTick = yTick;
axh1.YLim = yLim;
axh1.XTick = xTick;
axh1.XLim = xLim;
wantedTickLabel =  num2cell( yTick ) ;
axh1.YTickLabel =  wantedTickLabel;
shading(axh1,'interp');
xlabel(axh1,'x'); ylabel(axh1,'y') %rename x and y
axis square
% Polar order
axh2 = subplot(2,2,2); % Save the handle of the subplot
axh2.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh2);
h.TickLabelInterpreter = 'latex';
axh2.NextPlot = 'replaceChildren';
axh2.CLim = [0 1];
axh2.XLim = xLim; %row and columns are flipped
axh2.YLim = yLim; %row and columns are flipped
%axh2.YDir = 'rev';
axh2.XTick = xTick;
axh2.YTick = yTick;
axh2.YTickLabel =  wantedTickLabel;
shading(axh2,'interp');
xlabel(axh2,'x'); ylabel(axh2,'y') %rename x and y
axis square
% Distribution
axh3 = subplot(2,2,3); % Save the handle of the subplot
axh3.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh3);
h.Visible = 'off';
h.TickLabelInterpreter = 'latex';
maxDist = max(max( DistRec ) );
if maxDist <= 0
  maxDist = 0.1;
end
axh3.NextPlot = 'replaceChildren';
axh3.YLim = [0 maxDist];
axh3.XLim = [phi(1) phi(end)];
axis square
xlabel('$$\phi$$'); ylabel('f($$\phi$$)')
% Nematic order
axh4 = subplot(2,2,4); % Save the handle of the subplot
axh4.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh4);
h.TickLabelInterpreter = 'latex';
axh4.NextPlot = 'replaceChildren';
axh4.CLim = [0 1];
axh4.XLim = xLim; %row and columns are flipped
axh4.YLim = yLim; %row and columns are flipped
%axh4.YDir = 'rev';
axh4.XTick = xTick;
axh4.YTick = yTick;
axh4.YTickLabel = wantedTickLabel;
shading(axh4,'interp');
xlabel(axh4,'x'); ylabel(axh4,'y') %rename x and y
axis square
% Scale polar order by it's max value to for it changes.
polarTempX = OP.POPx_rec(SubInd1,SubInd2,:);
polarTempY = OP.POPy_rec(SubInd1,SubInd2,:);
maxPolar = max( max( max( OP.POP_rec(SubInd1,SubInd2,:) ) ) );
polarTempX = polarTempX ./ maxPolar;
polarTempY = polarTempY ./ maxPolar;
nemTempX = OP.NOPx_rec(SubInd1,SubInd2,:);
nemTempY = OP.NOPy_rec(SubInd1,SubInd2,:);
maxNem = max( max( max( OP.NOP_rec(SubInd1,SubInd2,:) ) ) );
nemTempX = nemTempX .* OP.NOP_rec(SubInd1,SubInd2,:) ./ maxNem;
nemTempY = nemTempY .* OP.NOP_rec(SubInd1,SubInd2,:) ./ maxNem;
% loop over frames
try
vec2loop = 1:nFrames;
  for ii = vec2loop
    % Concentration
    subplot(axh1);
    cla(axh1);
    pcolor( axh1, x, y, cScale * OP.C_rec(:,:,ii) );
    shading interp
    TitlStr = sprintf('Scale Concentration (bc) t = %.2f', TimeRec(ii));
    wantedTickLabel = flip( axh1.YTickLabel );
    axh1.YTickLabel =  wantedTickLabel;
    title(axh1,TitlStr);
    pause(0.001);
    drawnow;
    % Polar order
    subplot(axh2);
    cla(axh2);
    pcolor(axh2, x,y, OP.POP_rec(:,:,ii) );
    shading interp
    TitlStr = sprintf('Polar Order t = %.2f', TimeRec(ii));
    hold on
    quiver(axh2, x(SubInd2), y(SubInd1),...
      polarTempX(:,:,ii),...
      polarTempY(:,:,ii), 0,'color',[1,1,1] );
    title(axh2,TitlStr)
    pause(0.001);
    drawnow;
    % Distribution
    subplot(axh3);
    plot( axh3, phi, DistRec(:,ii) );
    TtlStr = sprintf('Distribution t = %.2f',...
      TimeRec(ii));
    title(axh3,TtlStr);
    pause(0.001);
    drawnow;
    % Nematic order
    cla(axh4)
    subplot(axh4);
    pcolor(axh4, x, y, OP.NOP_rec(:,:,ii) );
    shading interp
    TitlStr = sprintf('Nem. Order t = %.2f ', TimeRec(ii));
    hold on
    quiver(axh4, x(SubInd2),y(SubInd1),...
      nemTempX(:,:,ii), nemTempY(:,:,ii),0,...
      'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
    title(axh4,TitlStr)
    drawnow;
    pause(0.001);
    % Fig.Position can drift. Has to do with capturing a figure
    % before graphics render. Fix with drawnow and a pause.
    Fr = getframe(Fig);
    writeVideo(Mov,Fr);
  end% End frame loop
catch err
  fprintf('%s', err.getReport('extended')) ;
end % try catch
close(Fig)
