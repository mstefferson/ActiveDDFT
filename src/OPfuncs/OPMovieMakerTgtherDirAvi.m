
function OPMovieMakerTgtherDirAvi(MovStr,x,y,phi,OP,DistRec,TimeRec)
% Set up a indice vector so quiver is too crowded
Nx = length(x);
Ny = length(y);

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
FigWidth    = floor( ScreenHeight * .8 );
FigHeight   =  floor( ScreenHeight * .8);
FigPos      = [ floor( 0.5 * ( ScreenWidth - FigWidth ) ) ...
  floor( 0.5 * (ScreenHeight - FigHeight ) ) ...
  FigWidth FigHeight];

%Build a square box set by smallest dimension of screen

Fig = figure();
set(Fig, 'WindowStyle', 'normal');
Fig.Position = FigPos;

%Initialize the movie structure
Mov = VideoWriter(MovStr);
Mov.FrameRate = 4;
open(Mov);

%%%% Concentration %%%%%
nFrames = length(TimeRec);

set(gcf,'renderer','zbuffer')

axh1 = subplot(2,2,1); % Save the handle of the subplot
colorbar('peer',axh1)
axpos1 = get(axh1,'position'); % Save the position as ax
set(axh1,'NextPlot','replaceChildren',...
  'CLim',...
  [min(min(min(OP.C_rec))) max(max(max(OP.C_rec)))],...
  'YDir','normal');
axis square
set(axh1,'position',axpos1,'NextPlot','replaceChildren'); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')

axh2 = subplot(2,2,2); % Save the handle of the subplot
colorbar('peer',axh2)
axpos2 = get(axh2,'position'); % Save the position as ax
set(axh2,'position',axpos2,'NextPlot','replaceChildren',...
  'YDir','normal'); % Manually setting this holds the position with colorbar
axis square
xlabel('x'); ylabel('y')

% keyboard
axh3 = subplot(2,2,3); % Save the handle of the subplot
% colorbar('peer',axh3)
axpos3 = get(axh3,'position'); % Save the position as ax
set(axh3,'position',axpos3); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')
set(axh3,'NextPlot','replaceChildren','YLim',[ 0 max(max( DistRec ) )  ] ) ;
axis square
xlabel('\phi'); ylabel('f(\phi)')

axh4 = subplot(2,2,4); % Save the handle of the subplot
colorbar('peer',axh4)
axpos4 = get(axh4,'position'); % Save the position as ax
set(axh4,'position',axpos4,'NextPlot','replaceChildren'); % Manually setting this holds the position with colorbar
set(axh4,'NextPlot','replaceChildren',...
  'CLim',[0 1],'YDir','normal');
axis square
xlabel('x'); ylabel('y')

for ii = 1:nFrames
  try
  subplot(axh1);
  pcolor(axh1,x,y,OP.C_rec(:,:,ii)')
  shading(axh1,'interp');
  TitlStr = sprintf('Concentration t = %.2f', TimeRec(ii));
  title(axh1,TitlStr)
  
  % Polar order
  subplot(axh2);
  set(axh2,'NextPlot','replaceChildren',...
    'CLim',[0 max(max(max(OP.POP_rec)))],'YDir','normal');
  pcolor(axh2,x,y,OP.POP_rec(:,:,ii)');
  shading(axh2,'interp'); colorbar('peer',axh2);
  TitlStr = sprintf('Polar Order t = %.2f', TimeRec(ii));
  
  % What I and Matlab call x/y are switched
  hold on
  quiver(axh2,x(SubIndX),y(SubIndY),...
    OP.POPx_rec(SubIndX,SubIndY,ii)',...
    OP.POPy_rec(SubIndX,SubIndY,ii)' ,'color',[1,1,1]);
  hold off
  title(axh2,TitlStr)
   
  % Distribution
  subplot(axh3);
  plot( axh3, phi, DistRec(:,ii) );
  TtlStr = sprintf('Distribution t = %.2f',...
    TimeRec(ii));
  title(axh3,TtlStr);
  
  % Nematic order
  subplot(axh4);
  set(axh4,'NextPlot','replaceChildren',...
    'CLim',[0 1],'YDir','normal');
  pcolor(axh4, x, y, OP.NOP_rec(:,:,ii)');
  shading(axh4,'interp'); colorbar('peer',axh4)
  TitlStr = sprintf('Nem. Order t = %.2f ', TimeRec(ii));
  hold on
  quiver(axh4,x(SubIndX),y(SubIndY),...
    OP.NOPx_rec(SubIndX,SubIndY,ii)',OP.NOPy_rec(SubIndX,SubIndY,ii)',...
    'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
  hold off
  
  %     TitlStr = sprintf('Nematic Order');
  title(axh4,TitlStr)
  %     keyboard
  
  % Fig.Position can drift. Has to do with capturing a figure
  % before graphics render. Fix with drawnow and a pause
  drawnow;
  pause(0.001);
  
  Fr = getframe(Fig);
  writeVideo(Mov,Fr);
  catch
    keyboard
  end
  
end %% End frame loop


% keyboard
close all
