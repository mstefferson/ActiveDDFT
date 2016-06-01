function OPMovieMakerTgtherAvi(MovStr,x,y,phi,Crec,NOrec,POrec,DistRec,TimeRec)

% keyboard

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

% keyboard
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

axh2 = subplot(2,2,2); % Save the handle of the subplot
colorbar('peer',axh2)
axpos2 = get(axh2,'position'); % Save the position as ax
set(axh2,'position',axpos2); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')

% keyboard
axh3 = subplot(2,2,3); % Save the handle of the subplot
% colorbar('peer',axh3)
axpos3 = get(axh3,'position'); % Save the position as ax
set(axh3,'position',axpos3); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')

set(axh3,'NextPlot','replaceChildren','YLim',[ 0 max(max( DistRec ) )  ] ) ;

xlabel('\phi'); ylabel('f(\phi)')

axh4 = subplot(2,2,4); % Save the handle of the subplot
colorbar('peer',axh4)
axpos4 = get(axh4,'position'); % Save the position as ax
set(axh4,'position',axpos4); % Manually setting this holds the position with colorbar
xlabel('x'); ylabel('y')

% keyboard
for ii = 1:nFrames
%     keyboard
    %     set(axh1,'position',axpos1,'NextPlot','replaceChildren'); % Manually setting this holds the position with colorbar
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
    %     colorbar;
    %     set(gcf,'renderer','zbuffer')
    
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
    
    %     TitlStr = sprintf('Nematic Order');
    title(axh4,TitlStr)
    %     keyboard
    
%     Fr = getframe(Fig,[74 47 650 350]);
    Fr = getframe(Fig,[0 0 PosVec(3) PosVec(4)]);
    writeVideo(Mov,Fr);
    
%     if ii ~= nFrames
%         delete(axh1);
%         delete(axh2);
%         delete(axh3);
%         delete(axh4);
%     end
    
end %% End frame loop


% keyboard
close all
