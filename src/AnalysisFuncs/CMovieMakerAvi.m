% Makes movie of C vs time
function CMovieMakerAvi(MovStr,x,y,Crec,TimeRec)
% set font-size
fontSize = 34;
% Calculate log
if min( Crec(:) ) > 0
  logCrec = log( Crec );
  logLims = [min(min(min(logCrec))) max(max(max(logCrec)))];
else
  logCrec = ones( length(x), length(y), length(TimeRec) );
  logLims = [0 2];
end
% Set up figure
Fig = figure();
colormap(Fig, viridis);
set(Fig, 'WindowStyle', 'normal');
% Get screen size
screenSize = get(0,'screensize');
screenWidth = screenSize(3); screenHeight = screenSize(4);
figWidth    = floor( 0.86 * screenWidth );
figHeight   =  floor( 0.66 * screenHeight );
figPos      = [ floor( 0.5 * ( screenWidth - figWidth ) ) ...
  floor( 0.5 * (screenHeight - figHeight ) ) ...
  figWidth figHeight];
Fig.Position = figPos;
set(gcf,'renderer','zbuffer')
%Initialize the movie structure
Mov = VideoWriter(MovStr);
Mov.FrameRate = 4;
open(Mov);
%%%% Concentration %%%%%
nFrames = length(TimeRec);
subplot(1,2,1)
axh1 = gca; % Save the handle of the subplot
axis square
colorbar('peer',axh1)
axpos1 = get(axh1,'position'); % Save the position as ax
set(axh1,'NextPlot','replaceChildren',...
  'CLim', [min(min(min(Crec))) max(max(max(Crec)))],...
  'YDir','normal','Position',axpos1);
shading(axh1,'interp');
xlabel(axh1,'$$ x $$'); ylabel(axh1,'$$ y $$') 
axh1.FontSize = fontSize;
axis(axh1, 'square')
cTitle = '$$ C $$';
% log plot
subplot(1,2,2)
axh2 = gca; % Save the handle of the subplot
axis square
colorbar('peer',axh2)
axpos2 = get(axh2,'position'); % Save the position as ax
set(axh2,'NextPlot','replaceChildren',...
  'CLim', logLims ,...
  'YDir', 'normal','Position',axpos2);
shading(axh1,'interp');
xlabel(axh2,'$$ x $$'); ylabel(axh2,'$$ y $$') 
axh2.FontSize = fontSize;
axis(axh2, 'square')
clogTitle = '$$ \log{(C)} $$';
for ii = 1:nFrames
  % get t in title
  tTitle = [' ($$ t = $$ ' num2str( TimeRec(ii), '%.2f' ) ')' ];
  % Concentration
  subplot(axh1);
  pcolor(axh1,x,y,Crec(:,:,ii)')
  shading(axh1,'interp');
  TitlStr = [cTitle  tTitle]; 
  title(axh1,TitlStr)
  % Lot concentration
  subplot(axh2);
  pcolor(axh2,x,y,logCrec(:,:,ii)')
  shading(axh2,'interp');
  TitlStr = [clogTitle tTitle]; 
  title(axh2,TitlStr)
  % pause just in case
  drawnow 
  pause(0.1)
  % get frame and recprd
  Fr = getframe(Fig);
  %Fr = getframe(Fig);
  writeVideo(Mov,Fr);
end %% End frame loop
close all
