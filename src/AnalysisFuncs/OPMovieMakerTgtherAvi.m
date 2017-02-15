% Makes movie of OP vs time
function OPMovieMakerTgtherAvi(MovStr,x,y,phi,Crec,NOrec,POrec,DistRec,TimeRec)
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
axh1 = subplot(2,2,1); % Save the handle of the subplot
colorbar('peer',axh1)
axpos1 = get(axh1,'position'); % Save the position as ax
set(axh1,'NextPlot','replaceChildren',...
    'CLim',...
    [min(min(min(Crec))) max(max(max(Crec)))],...
    'YDir','normal');
set(axh1,'position',axpos1,'NextPlot','replaceChildren'); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')
% Polar Order
axh2 = subplot(2,2,2); % Save the handle of the subplot
colorbar('peer',axh2)
axpos2 = get(axh2,'position'); % Save the position as ax
set(axh2,'position',axpos2); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')
% Distribution
axh3 = subplot(2,2,3); % Save the handle of the subplot
axpos3 = get(axh3,'position'); % Save the position as ax
set(axh3,'position',axpos3); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')
set(axh3,'NextPlot','replaceChildren','YLim',[ 0 max(max( DistRec ) )  ] ) ;
xlabel('\phi'); ylabel('f(\phi)')
% Nematic
axh4 = subplot(2,2,4); % Save the handle of the subplot
colorbar('peer',axh4)
axpos4 = get(axh4,'position'); % Save the position as ax
set(axh4,'position',axpos4); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')
for ii = 1:nFrames
    % Concentration
    subplot(axh1);
    pcolor(axh1,x,y,Crec(:,:,ii)')
    shading(axh1,'interp');
    TitlStr = sprintf('Concentration t = %.2f', TimeRec(ii));
    title(axh1,TitlStr)
    % Polar order
    subplot(axh1);
    pcolor(axh2,x,y,POrec(:,:,ii)');
    shading(axh2,'interp'); colorbar('peer',axh2);
    set(axh2,'NextPlot','replaceChildren',...
        'CLim',[0 max(max(max(POrec)))],'YDir','normal');
    TitlStr = sprintf('Polar Order t = %.2f', TimeRec(ii));
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
    pcolor(axh4, x, y, NOrec(:,:,ii)');
    shading(axh4,'interp'); colorbar('peer',axh4)
    TitlStr = sprintf('Nem. Order t = %.2f ', TimeRec(ii));
    title(axh4,TitlStr)
    % get frame and recprd
    Fr = getframe(Fig,[0 0 PosVec(3) PosVec(4)]);
    writeVideo(Mov,Fr);
end %% End frame loop
close all
