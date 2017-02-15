% Makes movie of C vs time
function CMovieMakerAvi(MovStr,x,y,Crec,TimeRec)
% Set up figure
Fig = figure();
set(Fig, 'WindowStyle', 'normal');
PosVec = [680 558 1200 800];
Fig.Position = PosVec;
%Initialize the movie structure
Mov = VideoWriter(MovStr);
Mov.FrameRate = 4;
open(Mov);
%%%% Concentration %%%%%
nFrames = length(TimeRec);
set(gcf,'renderer','zbuffer')
axh1 = gca; % Save the handle of the subplot
colorbar('peer',axh1)
axpos1 = get(axh1,'position'); % Save the position as ax
set(axh1,'NextPlot','replaceChildren',...
  'CLim',...
  [min(min(min(Crec))) max(max(max(Crec)))],...
  'YDir','normal');
set(axh1,'position',axpos1,'NextPlot','replaceChildren'); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')
for ii = 1:nFrames
  % Concentration
  subplot(axh1);
  pcolor(axh1,x,y,Crec(:,:,ii)')
  shading(axh1,'interp');
  TitlStr = sprintf('Concentration t = %.2f', TimeRec(ii));
  title(axh1,TitlStr)
  % get frame and recprd
  Fr = getframe(Fig,[0 0 PosVec(3) PosVec(4)]);
  writeVideo(Mov,Fr);
end %% End frame loop
close all
