
function OPMovieMakerTgtherDirAvi(trial,x,y,phi,OP,DistRec,TimeRec)
% Set up a indice vector so quiver is too crowded
Nx = length(x);
Ny = length(y);

DivNumX = 8;
DivNumY = 10;
DeltaX  = floor(Nx / DivNumX );
DeltaY  = floor(Ny / DivNumY);
SubIndX = 1:DeltaX:(Nx + 1 - DeltaX);
SubIndY = 1:DeltaY:(Ny + 1 - DeltaY);

% Set up figure
Fig = figure();
set(Fig, 'WindowStyle', 'normal');
PosVec = [680 558 1200 800];
Fig.Position = PosVec;

%Initialize the movie structure
MovStr = sprintf('OPmovDir%.d.avi',trial);
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

% keyboard
for ii = 1:nFrames
%     keyboard
    %     set(axh1,'position',axpos1,'NextPlot','replaceChildren'); % Manually setting this holds the position with colorbarsubplot(1,2,2);
    subplot(axh1);
    pcolor(axh1,x,y,OP.C_rec(:,:,ii)')
    shading(axh1,'interp');
    TitlStr = sprintf('Concentration t = %.2f', TimeRec(ii));
    title(axh1,TitlStr)
    
    % Polar order
%     keyboard
    subplot(axh2);
    set(axh2,'NextPlot','replaceChildren',...
        'CLim',[0 max(max(max(OP.POP_rec)))],'YDir','normal');
    pcolor(axh2,x,y,OP.POP_rec(:,:,ii)');  
    shading(axh2,'interp'); colorbar('peer',axh2);
    TitlStr = sprintf('Polar Order t = %.2f', TimeRec(ii));
  
    % What I and Matlab call x/y are switched
    hold on
    quiver(axh2,x(SubIndX),y(SubIndY),...
        OP.nx_POP_rec(SubIndX,SubIndY,ii)',...
        OP.ny_POP_rec(SubIndX,SubIndY,ii)' ,'color',[1,1,1]);
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
        OP.NADx_rec(SubIndX,SubIndY,ii)',OP.NADy_rec(SubIndX,SubIndY,ii)',...
        'color',[1,1,1],'ShowArrowHead','off','LineWidth',0.1);
    hold off
   
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
%     keyboard
    
    
end %% End frame loop


% keyboard
close all
