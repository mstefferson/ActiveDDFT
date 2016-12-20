
function OPMovieMakerTgtherDirAvi(MovStr,x,y,phi,OP,DistRec,TimeRec)
nFrames = length(TimeRec);
Nx = length(x);
Ny = length(y);
Lx = x(end) + x(2);
Ly = y(end) + y(2);
% Find ticks
xTick = [0 Lx/2 Lx];
yTick = [0 Ly/2 Ly];
% Set up a index vector so quiver is too crowded
DivNumX = 8;
DivNumY = 8;
DeltaX  = ceil(Nx / DivNumX );
DeltaY  = ceil(Ny / DivNumY);
SubIndX = 1:DeltaX:(Nx + 1 - DeltaX);
SubIndY = 1:DeltaY:(Ny + 1 - DeltaY);
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
axh1.CLim = 1 /pi * [minC maxC];
axh1.XLim = [0 Lx]; %row and columns are flipped
axh1.YLim = [0 Ly]; %row and columns are flipped
axh1.YDir = 'rev';
% axh1.YDir = 'normal';
axh1.XTick = xTick;
axh1.YTick = yTick;
% wantedTickLabel =  axh1.YTickLabel;
wantedTickLabel = flip( axh1.YTickLabel );
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
axh2.XLim = [0 Lx]; %row and columns are flipped
axh2.YLim = [0 Ly]; %row and columns are flipped
axh2.YDir = 'rev';
% axh2.YDir = 'normal';
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
axh3.XLim = [0 2*pi];
axis square
xlabel('$$\phi$$'); ylabel('f($$\phi$$)')
% Nematic order
axh4 = subplot(2,2,4); % Save the handle of the subplot
axh4.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh4);
h.TickLabelInterpreter = 'latex';
axh4.NextPlot = 'replaceChildren';
axh4.CLim = [0 1];
axh4.XLim = [0 Lx]; %row and columns are flipped
axh4.YLim = [0 Ly]; %row and columns are flipped
% axh4.YDir = 'normal';
axh4.YDir = 'rev';
axh4.XTick = xTick;
axh4.YTick = yTick;
axh4.YTickLabel = wantedTickLabel;
shading(axh4,'interp');
xlabel(axh4,'x'); ylabel(axh4,'y') %rename x and y
axis square
% Scale polar order by it's max value to for it changes.
polarTempX = OP.POPx_rec(SubIndX,SubIndY,:);
polarTempY = OP.POPy_rec(SubIndX,SubIndY,:);
maxPolar = max( max( max( OP.POP_rec(SubIndX,SubIndY,:) ) ) );
polarTempX = polarTempX ./ maxPolar;
polarTempY = polarTempY ./ maxPolar;
nemTempX = OP.NOPx_rec(SubIndX,SubIndY,:);
nemTempY = OP.NOPy_rec(SubIndX,SubIndY,:);
maxNem = max( max( max( OP.NOP_rec(SubIndX,SubIndY,:) ) ) );
nemTempX = nemTempX .* OP.NOP_rec(SubIndX,SubIndY,:) ./ maxNem;
nemTempY = nemTempY .* OP.NOP_rec(SubIndX,SubIndY,:) ./ maxNem;
% keyboard
% loop over frames
try
 vec2loop = nFrames;
%  vec2loop = 1:nFrames;
  for ii = vec2loop
    % Concentration
    subplot(axh1);
    cla(axh1);
    imagesc(axh1, x, y, rot90( OP.C_rec(:,:,ii) ) ./ pi);
    TitlStr = sprintf('Scale Concentration (bc) t = %.2f', TimeRec(ii));
    title(axh1,TitlStr);
    pause(0.001);
    drawnow;
    % Polar order
    subplot(axh2);
    cla(axh2);
    imagesc(axh2,x,y, rot90( OP.POP_rec(:,:,ii) ) );
    TitlStr = sprintf('Polar Order t = %.2f', TimeRec(ii));
    hold on
    quiver(axh2,x(SubIndX),y(SubIndY),...
      rot90( polarTempX(:,:,ii) ),...
      rot90( -polarTempY(:,:,ii) ), 0,'color',[1,1,1]);
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
    imagesc(axh4, x, y, rot90( OP.NOP_rec(:,:,ii) ) );
    TitlStr = sprintf('Nem. Order t = %.2f ', TimeRec(ii));
    hold on
    quiver(axh4, x(SubIndX),y(SubIndY),...
      rot90( nemTempX(:,:,ii) ), rot90( -nemTempY(:,:,ii) ),0,...
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
