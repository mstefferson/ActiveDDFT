
function OPMovieMakerTgtherDirAvi(MovStr,x,y,phi,OP,DistRec,TimeRec)
% Set up a indice vector so quiver is too crowded
Nx = length(x);
Ny = length(y);
% plot axis divisions
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
% set(Fig, 'WindowStyle', 'normal');
Fig.Position = FigPos;
%Initialize the movie structure
Mov = VideoWriter(MovStr);
Mov.FrameRate = 4;
open(Mov);
% C
nFrames = length(TimeRec);
set(gcf,'renderer','zbuffer')
axh1 = subplot(2,2,1); % Save the handle of the subplot
axh1.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh1);
h.TickLabelInterpreter = 'latex';
minC = min(min(min(OP.C_rec )));
maxC = max(max(max(OP.C_rec)));
if minC >= maxC
  minC = 0;
  maxC = 0.1;
end
axh1.NextPlot = 'replaceChildren';
axh1.CLim = 1 /pi * [minC maxC];
% axh1.YDir = 'normal';
axis square
% xlabel('x'); ylabel('y')
% Polar order
axh2 = subplot(2,2,2); % Save the handle of the subplot
axh2.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh2);
h.TickLabelInterpreter = 'latex';
axh2.NextPlot = 'replaceChildren';
axh2.CLim = [0 1];
% axh2.YDir = 'normal';
axis square
% xlabel('x'); ylabel('y')
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
axis square
xlabel('$$\phi$$'); ylabel('f($$\phi$$)')
% Nematic order
axh4 = subplot(2,2,4); % Save the handle of the subplot
axh4.TickLabelInterpreter = 'latex';
h = colorbar('peer',axh4);
h.TickLabelInterpreter = 'latex';
axh4.NextPlot = 'replaceChildren';
axh4.CLim = [0 1];
% axh4.YDir = 'normal';
axh4.YTick = linspace(0, Nx, 10)  ;

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
% loop over framces
% keyboard
try

%   for ii = 1:nFrames
    ii = nFrames;
    % Concentration
    subplot(axh1);
%     pcolor(axh1,x,y,OP.C_rec(:,:,ii)' ./ pi)
%     imagesc(axh1,y,x,OP.C_rec(:,:,ii) ./ pi)
    % Matlab's way
    imagesc(axh1, y, x, OP.C_rec(:,:,ii) ./ pi)
    axh1.YDir = 'rev';

    shading(axh1,'interp');
    TitlStr = sprintf('Scale Concentration (bc) t = %.2f', TimeRec(ii));
    title(axh1,TitlStr)
    drawnow;
    pause(0.001);
    % Polar order
    subplot(axh2);
    set(axh2,'NextPlot','replaceChildren',...
      'CLim',[0 1]);
    axh2.YDir = 'rev';
    %pcolor(axh2,x,y,OP.POP_rec(:,:,ii)');
    %imagesc(axh2,y,x,OP.POP_rec(:,:,ii));
    %imageQuiver( amps, xDir, yDir, myMethod, ax)
    shading(axh2,'interp'); 
    TitlStr = sprintf('Polar Order t = %.2f', TimeRec(ii));
    % What I and Matlab call x/y are switched
    % matlab's way
    %imagesc(axh2,OP.POP_rec(:,:,ii));
    %quiver(axh2,y(SubIndX),x(SubIndY),...
      %polarTempX(:,:,ii),...
      %polarTempY(:,:,ii), 0,'color',[1,1,1]);
    hold on
    quiver(axh2,y(SubIndY),x(SubIndX),...
      polarTempX(:,:,ii),...
      -polarTempY(:,:,ii), 0,'color',[1,1,1]);
    hold off
    %imageQuiver( OP.POP_rec(:,:,ii), OP.POPx_rec(:,:,ii), ...
    %  OP.POPy_rec(:,:,ii),1, axh2 )

    title(axh2,TitlStr)
    drawnow;
    pause(0.001);
    
    % Distribution
    subplot(axh3);
    plot( axh3, phi, DistRec(:,ii) );
    TtlStr = sprintf('Distribution t = %.2f',...
      TimeRec(ii));
    title(axh3,TtlStr);
    drawnow;
    pause(0.001);
    % Nematic order
    subplot(axh4);
    set(axh4,'NextPlot','replaceChildren',...
      'CLim',[0 1]);
    %pcolor(axh4, x, y, OP.NOP_rec(:,:,ii)');
    % Matlab's way
    imageQuiver( OP.POP_rec(:,:,ii), OP.POPx_rec(:,:,ii), ...
      OP.POPy_rec(:,:,ii),1, axh2 )
    %imagesc(axh4, y, x, OP.NOP_rec(:,:,ii));
    imagesc(axh4,y, x, OP.NOP_rec(:,:,ii));
    axh4.YDir = 'rev';
    shading(axh4,'interp'); colorbar('peer',axh4)
    TitlStr = sprintf('Nem. Order t = %.2f ', TimeRec(ii));
    hold on
    %quiver(axh4,y(SubIndX),x(SubIndY),...
      %nemTempX(:,:,ii)',nemTempY(:,:,ii)',0,...
      %'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
    quiver(axh4,y(SubIndY),x(SubIndX),...
      nemTempX(:,:,ii),-nemTempY(:,:,ii)',0,...
      'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
    hold off
    title(axh4,TitlStr)
    drawnow;
    pause(0.001);
    % Fig.Position can drift. Has to do with capturing a figure
    % before graphics render. Fix with drawnow and a pause.
    Fr = getframe(Fig);
    writeVideo(Mov,Fr);
    keyboard
%   end% End frame loop
  
catch err
  fprintf('%s', err.getReport('extended')) ;
  %   keyboard
end % try catch

close(Fig)
