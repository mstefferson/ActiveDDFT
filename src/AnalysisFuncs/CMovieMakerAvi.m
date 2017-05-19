% Makes movie of C vs time
function CMovieMakerAvi(MovStr,x,y,Crec,TimeRec)
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
set(Fig, 'WindowStyle', 'normal');
% Get screen size
pixSS = get(0,'screensize');
pixW = pixSS(3);
pixH = pixSS(4);
figSize = [pixW pixH] / 2;
figCorner = [ ( pixW-figSize(1) ) ( pixH-figSize(2) ) ] ./ 2;
posVec = [figCorner figSize];
Fig.Position = posVec;
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
xlabel('x'); ylabel('y')
% log plot
subplot(1,2,2)
axh2 = gca; % Save the handle of the subplot
axis square
colorbar('peer',axh2)
axpos2 = get(axh2,'position'); % Save the position as ax
set(axh2,'NextPlot','replaceChildren',...
  'CLim', logLims ,...
  'YDir', 'normal','Position',axpos2);
xlabel('x'); ylabel('y');
for ii = 1:nFrames
  % Concentration
  subplot(axh1);
  pcolor(axh1,x,y,Crec(:,:,ii)')
  shading(axh1,'interp');
  TitlStr = sprintf('Concentration t = %.2f', TimeRec(ii));
  title(axh1,TitlStr)
  % Lot concentration
  subplot(axh2);
  pcolor(axh2,x,y,logCrec(:,:,ii)')
  shading(axh2,'interp');
  TitlStr = sprintf('Log Concentration t = %.2f', TimeRec(ii));
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
