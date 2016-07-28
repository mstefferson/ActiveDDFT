Nm = PTMGDObj.ParamObj.Nm;
AngMode = Nm/2;
nFrames = length(DenRecObj.TimeRecVec);
dt_rec  = PTMGDObj.TimeObj.t_record;
dPhi    = 2 / Nm * 180;

movie      = 1;
DenAll     = 0;
DenPosFixed   = 1;
DenLast    = 0;
InterAll   = 0;
InterFixed = 0;
InterDenFixed = 0;
InterLast  = 0;


%%%%%%%%%%%%%%%% Watch the density of the individual modes %%%%%%%%%%%%%

if DenAll
    Frames2RunEnd  = nFrames;
    Frames2RunStrt = 1;
    %     Frames2Run = 5;
    
    figure()
    set(gca,'NextPlot','replaceChildren',...
        'CLim',[ min(min(min(min( DenRecObj.Density_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))...
        max(max(max(max( DenRecObj.Density_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) )))) ],...
        'ZLim',[ min(min(min(min( DenRecObj.Density_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))...
        max(max(max(max( DenRecObj.Density_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) )))) ],...
        'YDir','normal');
    set(gcf,'renderer','zbuffer')
    
    
    
    for ii = 1:Nm
        for jj = Frames2RunStrt : Frames2RunEnd
            pause(.2)
            TtlStr = sprintf('Phi = %f t = %f', dPhi*(ii-1), dt_rec*(jj-1));
            surf(DenRecObj.Density_rec(:,:,ii,jj));
            title(TtlStr);
            getframe(gcf);
        end
    end
end

% See what the angular distribution looks like at a fixed gridpoint
if DenPosFixed
    figure()
    xPos  = Nm/2;
    yPos  = Nm/2;
    % Frames2Run = nFrames;
    Frames2RunEnd  = nFrames;
    Frames2RunStrt = 1;
    %
    % Frames2RunEnd  = 170;
    % Frames2RunStrt = 140;
    set( gca,'NextPlot','replaceChildren',...
        'YLim',[ min(min(min(min( DenRecObj.Density_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ...
        max(max(max(max( DenRecObj.Density_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ] )
    set(gcf,'renderer','zbuffer')
    if movie
       MovM_BT(1:( Frames2RunEnd + 1 - Frames2RunStrt) ) = ...
           struct('cdata',[], 'colormap',[]);
    end
    
    for ii = Frames2RunStrt  : Frames2RunEnd
        pause(.1);
        TtlStr = sprintf('t = %f', dt_rec*(ii-1));
        plot( PTMGDObj.GridObj.phi, reshape( DenRecObj.Density_rec(xPos,yPos,:,ii), 1, Nm ) );
        title(TtlStr);
        if movie
%             keyboard
            MovM_BT(ii) = getframe(gcf);
        else
            getframe(gcf)
        end
    end

end

% Cycle through the angles on the last frame

if DenLast
    figure()
    set( gca,'NextPlot','replaceChildren',...
        'YLim',[ min(min(min(min( DenRecObj.Density_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ...
        max(max(max(max( DenRecObj.Density_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ] )
    
    set(gcf,'renderer','zbuffer')
    for ii = 1:Nm
        TtlStr = sprintf('Phi = %f t = %f', dPhi*(ii-1), dt_rec*(nFrames-1));
        surf(DenRecObj.Density_rec(:,:,ii,nFrames));
        title(TtlStr);
        getframe(gcf);
    end
    
end
%%%%%%%%%%%%%%%% Watch the dRho from interactions of the individual modes %%%%%%%%%%%%%
if InterAll
    Frames2RunEnd  = nFrames;
    Frames2RunStrt = 1;
    %
    % Frames2RunEnd  = 170;
    % Frames2RunStrt = 160;
    
    set(gca,'NextPlot','replaceChildren',...
        'CLim',[min(min(min(min( DenRecObj.dRhoInter_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))...
        max(max(max(max( DenRecObj.dRhoInter_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))],...
        'ZLim',[min(min(min(min( DenRecObj.dRhoInter_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))...
        max(max(max(max( DenRecObj.dRhoInter_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))],...
        'YDir','normal');
    set(gcf,'renderer','zbuffer')
    
    for ii = 1:Nm
        for jj = Frames2RunStrt  : Frames2RunEnd
            pause(.01)
            TtlStr = sprintf('Phi = %f t = %f', dPhi*(ii-1), dt_rec*(jj-1));
            surf(DenRecObj.dRhoInter_rec(:,:,ii,jj));
            title(TtlStr);
            getframe(gcf);
        end
    end
    
end
% See what drho from interactions looks like at a fixed gridpoint
if InterFixed
    figure()
    xPos  = Nm/2;
    yPos  = Nm/2;
    % Frames2Run = nFrames;
    Frames2RunEnd  = nFrames;
    Frames2RunStrt = 1;
    
    % Frames2RunEnd  = 8;
    % Frames2RunStrt = 130;
    
    set( gca,'NextPlot','replaceChildren',...
        'YLim',[ min(min(min(min( DenRecObj.dRhoInter_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ...
        max(max(max(max( DenRecObj.dRhoInter_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ] )
    
    set(gcf,'renderer','zbuffer')
    
    for ii = Frames2RunStrt  : Frames2RunEnd
        pause(.1);
        TtlStr = sprintf('t = %f', dt_rec*(ii-1));
        plot( ParaTimMemObj.GridObj.phi, reshape( DenRecObj.dRhoInter_rec(xPos,yPos,:,ii), 1, Nm ) );
        title(TtlStr);
        getframe(gcf);
    end
    
    dRhoInt_FT = fftshift(fftn(DenRecObj.dRhoInter_rec(:,:,:,2)));
end

% See what drho from interactions and density looks like at a fixed gridpoint
if InterDenFixed
    
    xPos  = Nm/2;
    yPos  = Nm/2;
    Frames2RunEnd  = nFrames;
    Frames2RunStrt = 1;
    
    figure()
    set(gcf,'renderer','zbuffer')
    subplot(2,1,1)
    set( gca,'NextPlot','replaceChildren',...
        'YLim',[ min(min(min(min( DenRecObj.Density_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ...
        max(max(max(max( DenRecObj.Density_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ] )
    
    subplot(2,1,2)
    set( gca,'NextPlot','replaceChildren',...
        'YLim',[ min(min(min(min( DenRecObj.dRhoInter_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ...
        max(max(max(max( DenRecObj.dRhoInter_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ] )
    
    for ii = Frames2RunStrt  : Frames2RunEnd
        %         pause(.1);
        TtlStr = sprintf('Density t = %f', dt_rec*(ii-1));
        subplot(2,1,1)
        %         plot( ParaTimMemObj.GridObj.phi, reshape( DenRecObj.Density_rec(xPos,yPos,:,ii), 1, N ) );
        plot( reshape( DenRecObj.Density_rec(xPos,yPos,:,ii), 1, Nm ) );
        title(TtlStr);
        subplot(2,1,2)
        TtlStr = sprintf('dRho t = %f', dt_rec*(ii-1));
        %         plot( ParaTimMemObj.GridObj.phi, reshape( DenRecObj.dRhoInter_rec(xPos,yPos,:,ii), 1, N ) );
        plot( reshape( DenRecObj.dRhoInter_rec(xPos,yPos,:,ii), 1, Nm ) );
        title(TtlStr);
        getframe(gcf);
    end
    
end

% Cycle through the angles on the last frame

if InterLast
    for ii = 1:Nm
        TtlStr = sprintf('Phi = %f t = %f', dPhi*(ii-1), dt_rec*(nFrames-1));
        surf(DenRecObj.dRhoInter_rec(:,:,ii,nFrames));
        title(TtlStr);
        getframe(gcf);
    end
end
% Cycle through the angles on the last frame


if 0
    %%%%%%%%%%%%%%%% Watch the dRho from propagator of the individual modes %%%%%%%%%%%%%
    set(gca,'NextPlot','replaceChildren',...
        'CLim',[min(min(min(min( DenRecObj.dRhoProp_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))...
        max(max(max(max( DenRecObj.dRhoProp_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))],...
        'ZLim',[min(min(min(min( DenRecObj.dRhoProp_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))...
        max(max(max(max( DenRecObj.dRhoProp_rec( :, :, :, Frames2RunStrt:Frames2RunEnd ) ))))],...
        'YDir','normal');
    set(gcf,'renderer','zbuffer')
    
    for ii = 1:N
        for jj = 1:nFrames
            TtlStr = sprintf('Phi = %f t = %f', dPhi*(ii-1), dt_rec*(jj-1));
            surf(DenRecObj.dRhoProp_rec(:,:,ii,jj));
            title(TtlStr);
            getframe(gcf);
        end
    end
    
    % See what drho from propagator looks like at a fixed gridpoint
    figure()
    xPos  = N/2;
    yPos  = N/2;
    % Frames2Run = nFrames;
    Frames2RunEnd  = nFrames;
    Frames2RunStrt = 1;
    
    % Frames2RunEnd  = 170;
    % Frames2RunStrt = 130;
    set( gca,'NextPlot','replaceChildren',...
        'YLim',[ min(min(min(min( DenRecObj.dRhoProp_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ...
        max(max(max(max( DenRecObj.dRhoProp_rec( xPos, yPos, :, Frames2RunStrt:Frames2RunEnd ) )))) ] )
    
    set(gcf,'renderer','zbuffer')
    
    
    for ii = Frames2RunStrt  : Frames2RunEnd
        %     pause(.1);
        TtlStr = sprintf('t = %f', dt_rec*(ii-1));
        plot( ParaTimMemObj.GridObj.phi, reshape( DenRecObj.dRhoProp_rec(xPos,yPos,:,ii), 1, N ) );
        title(TtlStr);
        getframe(gcf);
    end
    
    % Cycle through the angles on the last frame
    
    for ii = 1:N
        TtlStr = sprintf('Phi = %f t = %f', dPhi*(ii-1), dt_rec*(nFrames-1));
        surf(DenRecObj.dRhoProp_rec(:,:,ii,nFrames));
        title(TtlStr);
        getframe(gcf);
    end
    
    Mode = 2;
    Frame = 1;
    figure
    subplot(3,1,1)
    surf(DenRecObj.Density_rec(:,:,Mode,Frame))
    title('Density')
    
    subplot(3,1,2)
    surf(DenRecObj.dRhoInter_rec(:,:,Mode,Frame))
    title('Density change from interactions')
    
    subplot(3,1,3)
    surf(DenRecObj.dRhoProp_rec(:,:,Mode,Frame))
    title('Density change from propagator')
    
    
end




