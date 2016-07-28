function [MovieObj] = ...
    OPMatMovieMakerSep(Nx,Ny,Nm,x,y,phi,OrderParamObj,Density_rec,trial)

%%%% Concentration %%%%%

nFrames = OrderParamObj.nFrames;

%Initialize the movie structure array
M_C(nFrames) = struct('cdata',zeros(Nx,Ny,3,'int8'), 'colormap',[]); %initialize movie stucture
% Set up figure

set(gca,'NextPlot','replaceChildren','CLim',[min(min(min(OrderParamObj.C_rec))) max(max(max(OrderParamObj.C_rec)))],'YDir','normal');
set(gcf,'renderer','zbuffer')
% set(gcf,'renderer','zbuffer')
colorbar

for ii = 1:nFrames
    pcolor(x,y,OrderParamObj.C_rec(:,:,ii)')
    shading interp;
    TitlStr = sprintf('Concentration t = %f', OrderParamObj.TimeRec(ii));
    title(TitlStr)
    %     keyboard
    %change the scale
    M_C(ii) = getframe(gcf); %Store the frame
end

close all

% keyboard
%%%%Polar order%%%%%%

% Initialize movie stucture array
M_P(nFrames) = struct('cdata',zeros(Nx,Ny,3,'int8'), 'colormap',[]);

%Set up figure
figure()
xlabel('x')
ylabel('y')
set(gca,'NextPlot','replaceChildren','CLim',[0 max(max(max(OrderParamObj.POP_rec)))],'YDir','normal');
set(gcf,'renderer','zbuffer')
for ii = 1:nFrames
    set(gca,'NextPlot','replaceChildren'); %Because of hold, need to do this again
    pcolor(x,y,OrderParamObj.POP_rec(:,:,ii)');
    shading interp;
    colorbar;
    %     set(gcf,'renderer','zbuffer')
    hold on
    quiver(x,y,OrderParamObj.nx_POP_rec(:,:,ii)',OrderParamObj.ny_POP_rec(:,:,ii)','color',[1,1,1]);
    hold off
    TitlStr = sprintf('Polar Order t = %f', OrderParamObj.TimeRec(ii));
    title(TitlStr)
    M_P(ii) = getframe(gcf); %Store the frame
    %         keyboard
end

close all
% keyboard
%%%% Nematic order%%%%%%

% Initialize movie stucture array
M_N(nFrames) = struct('cdata',zeros(Nx,Ny,3,'int8'), 'colormap',[]);

% Set up figure
xlabel('x')
ylabel('y')
set(gca,'NextPlot','replaceChildren','CLim',[0 max(max(max(OrderParamObj.NOP_rec)))],'YDir','normal');
set(gcf,'renderer','zbuffer')
for ii = 1:nFrames
    set(gca,'NextPlot','replaceChildren'); %Because of hold, need to do this again
    pcolor(x,y,OrderParamObj.NOP_rec(:,:,ii)');
    shading interp;
    colorbar;
    hold on
    quiver(x,y,OrderParamObj.NADx_rec(:,:,ii)',OrderParamObj.NADy_rec(:,:,ii)',...
        'color',[1,1,1],'ShowArrowHead','off');
    hold off
    TitlStr = sprintf('Nematic Order t = %f', OrderParamObj.TimeRec(ii));
    title(TitlStr)
    M_N(ii) = getframe(gcf); %Store the frame
    %         keyboard
end

%%%%%Distribution at fized position%%%%%%%%
xPos  = Nx/2;
yPos  = Ny/2;
%
% Frames2RunEnd  = 170;
% Frames2RunStrt = 140;
% keyboard
eps = 0.000001;
set( gca,'NextPlot','replaceChildren',...
    'YLim',[ min(min(min(min( Density_rec( xPos, yPos, :, : ) ) ) ) ) - eps ...
    max(max(max(max( Density_rec( xPos, yPos, :, : ) ) ) ) )  + eps] )
set(gcf,'renderer','zbuffer')

M_BT(1:nFrames ) = ...
    struct('cdata',[], 'colormap',[]);

for ii = 1  : nFrames
    
    TtlStr = sprintf('t = %f xind = %d yind = %d ', OrderParamObj.TimeRec(ii), xPos,yPos);
    plot( phi, reshape( Density_rec(xPos,yPos,:,ii), 1, Nm ) );
    title(TtlStr);
    M_BT(ii) = getframe(gcf);
    
end
% keyboard
MovieObj = struct('M_C',M_C,'M_P',M_P,'M_N',M_N,'M_BT',M_BT);


        
        MovM_C  = MovieObj.M_C;
        MovM_P  = MovieObj.M_P;
        MovM_N  = MovieObj.M_N;
        MovM_BT = MovieObj.M_BT;
        
        
        MovMCStr = sprintf('MovieMC_%i',trial);
        MovMPStr = sprintf('MovieMP_%i',trial);
        MovMNStr = sprintf('MovieMN_%i',trial);
        MovBTStr = sprintf('MovieBT_%i',trial);
        save(MovMCStr,'MovM_C','-v7.3')
        save(MovMNStr,'MovM_N','-v7.3')
        save(MovMPStr,'MovM_P','-v7.3')
        save(MovBTStr,'MovM_BT','-v7.3')
        
close all