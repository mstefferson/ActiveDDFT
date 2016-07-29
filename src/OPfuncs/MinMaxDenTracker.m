% Function: MinMaxDenTracker.m
% Plots the min/max density at each spatial gridpoint
%
% Input:
% MinMaxDenTracker(DenRecObj.Density_rec,DenRecObj.TimeRecVec,PTMGDObj.ParamObj,PTMGDObj.gridObj)
% Output:
% A movie!!!

function [M_MinMax,M_MinMaxSS,M_dist] = MinMaxDenTracker(Density_rec,TimeRecVec,ParamObj,gridObj)

MinMaxPlot         = 1;
MinMaxSSPlot       = 1;
AdjacentPointsPlot = 0;

%Initialize some things
nFrames    = length(TimeRecVec) - 1;
M_dist(nFrames) = ...
    struct('cdata',[], 'colormap',[]);
M_MinMax(nFrames) = ...
    struct('cdata',zeros(systemObj.Nx,systemObj.Ny,3,'int8'), 'colormap',[]); %initialize movie stucture
M_MinMaxSS(nFrames) = struct('cdata',zeros(systemObj.Nx,systemObj.Ny,3,'int8'), ...
    'colormap',[]); %initialize movie stucture

if MinMaxPlot
    phiMaxSpat = zeros(systemObj.Nx,systemObj.Ny);
    phiMinSpat = phiMaxSpat;
    AbsMinRho = min( min( min( min( Density_rec ) ) ) );
    AbsMaxRho = max( max( max( max( Density_rec ) ) ) );
    
    % Movie stuff
    
    % Set up figure
    
    figure(1)
    set(gca,'NextPlot','replaceChildren',...
        'CLim',[ min(min(min(min( Density_rec )))) max(max(max(max( Density_rec )))) ],...
        'YDir','normal');
    set(gcf,'renderer','zbuffer')
    % set(gcf,'renderer','zbuffer')
    colorbar
    
    % % keyboard
    for t = 1:nFrames
        
        rhoMax    = max( Density_rec(:,:,:,t),[],3);
        rhoMin    = min( Density_rec(:,:,:,t),[],3);
        NumPart = trapz_periodic( trapz_periodic( trapz_periodic( Density_rec(:,:,:,t) ) ) );
        for i = 1 : systemObj.Nx
            for j = 1 : systemObj.Ny
                %             keyboard
                % Find max phi
                
                [~, IndMaxPhiTempFixedSpat] =  max( Density_rec(i,j,:,t) );
                if length( IndMaxPhiTempFixedSpat ) > 1
                    IndMaxPhiTempFixedSpat = IndMaxPhiTempFixedSpat(1);
                end
                % Find min phi
                [~, IndMinPhiTempFixedSpat] =  min( Density_rec(i,j,:,t) );
                if length( IndMinPhiTempFixedSpat ) > 1
                    IndMinPhiTempFixedSpat = IndMinPhiTempFixedSpat(1);
                end
                
                % Easy if we just keep angles between [0, pi/2]
                % We can be more careful about this later
                if  IndMaxPhiTempFixedSpat > systemObj.Nm / 2
                    %                 keyboard
                    IndMaxPhiTempFixedSpat = IndMaxPhiTempFixedSpat - systemObj.Nm / 2;
                end
                
                
                if  IndMinPhiTempFixedSpat > systemObj.Nm / 2
                    IndMinPhiTempFixedSpat = IndMinPhiTempFixedSpat - systemObj.Nm / 2;
                end
                
                phiMaxSpat(i,j)     = gridObj.phi( IndMaxPhiTempFixedSpat );
                phiMinSpat(i,j)     = gridObj.phi( IndMinPhiTempFixedSpat );
            end
        end
        
        %      keyboard
        subplot(3,1,1)
        pcolor(gridObj.x,gridObj.y,rhoMax')
        set(gca,'CLim', [AbsMinRho AbsMaxRho],'YDir','rev');
        shading interp;
        colorbar;
        hold on
        quiver(gridObj.x,gridObj.y,cos(phiMaxSpat), sin(phiMaxSpat),...
            'color',[0 0 0],'AutoScaleFactor',0.5 );
        hold off
        TitlStr = sprintf('Min Max Plot t = %f N~%f', TimeRecVec(t), NumPart );
        title(TitlStr)
        
        subplot(3,1,2)
        pcolor(gridObj.x,gridObj.y,rhoMin');
        set(gca,'CLim', [AbsMinRho AbsMaxRho],'YDir','rev');
        shading interp;
        colorbar;
        hold on
        quiver(gridObj.x,gridObj.y,cos(phiMinSpat)', sin(phiMinSpat)',...
            'color',[0 0 0],'AutoScaleFactor',0.5 );
        hold off
        
        subplot(3,1,3)
        pcolor(gridObj.x,gridObj.y, mean( Density_rec(:,:,:,t),3 )' );
        set(gca,'CLim', [AbsMinRho AbsMaxRho],'YDir','rev');
        shading interp;
        colorbar;
        
        % Store Frame
        M_MinMax(t) =  getframe(gcf);
    end %end frame loop
    
end %MinMaxPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MinMaxSSPlot
    % Just Plot a subset of the
    NxSs = 8:24;
    NySs = 8:24;
    M_MinMaxSS(nFrames) = struct('cdata',zeros( length(NxSs), length(NySs) ,3,'int8'), ...
        'colormap',[]); %initialize movie stucture
    
    % keyboard
    phiMaxSpat = zeros( length(NxSs), length(NySs) );
    phiMinSpat = phiMaxSpat;
    AbsMinRho = min( min( min( min( Density_rec(NxSs,NySs,:,:) ) ) ) );
    AbsMaxRho = max( max( max( max( Density_rec(NxSs,NySs,:,:) ) ) ) );
    
    % Set up figure
    
    figure(2)
    set(gca,'NextPlot','replaceChildren',...
        'CLim',[ AbsMinRho AbsMaxRho ],...
        'YDir','normal');
    set(gcf,'renderer','zbuffer')
    % set(gcf,'renderer','zbuffer')
    colorbar
    
    % keyboard
    for t = 1:nFrames
        
        rhoMax    = max( Density_rec(NxSs,NySs,:,t), [],3);
        rhoMin    = min( Density_rec(NxSs,NySs,:,t), [],3);
        NumPart = trapz_periodic( trapz_periodic( trapz_periodic( Density_rec(NxSs,NySs,:,t) ) ) );
        for i = 1:length(NxSs)
            for j = 1:length(NySs)
                %             keyboard
                
                % Find max phi
                [~, IndMaxPhiTempFixedSpat] =  max( Density_rec(NxSs(i),NySs(j),:,t) );
                if length( IndMaxPhiTempFixedSpat ) > 1
                    IndMaxPhiTempFixedSpat = IndMaxPhiTempFixedSpat(1);
                end
                % Find min phi
                [~, IndMinPhiTempFixedSpat] =  min( Density_rec(NxSs(i),NySs(j),:,t) );
                
                if length( IndMinPhiTempFixedSpat ) > 1
                    IndMinPhiTempFixedSpat = IndMinPhiTempFixedSpat(1);
                end
                
                % Easy if we just keep angles between [0, pi/2]
                % We can be more careful about this later
                if  IndMaxPhiTempFixedSpat > systemObj.Nm / 2
                    IndMaxPhiTempFixedSpat = IndMaxPhiTempFixedSpat - systemObj.Nm / 2;
                end
                
                
                if  IndMinPhiTempFixedSpat > systemObj.Nm / 2
                    IndMinPhiTempFixedSpat = IndMinPhiTempFixedSpat - systemObj.Nm / 2;
                end
                
                phiMaxSpat(i,j)     = gridObj.phi( IndMaxPhiTempFixedSpat );
                phiMinSpat(i,j)     = gridObj.phi( IndMinPhiTempFixedSpat );
                
                gridObj.phi( IndMaxPhiTempFixedSpat );
                gridObj.phi( IndMinPhiTempFixedSpat );
                %             keyboard
                
            end
        end
        
        %      keyboard
        
        subplot(3,1,1)
        pcolor(gridObj.x(NxSs),gridObj.y(NySs),rhoMax')
        set(gca,'CLim', [AbsMinRho AbsMaxRho],'YDir','rev');
        shading interp;
        colorbar;
        hold on
        quiver(gridObj.x(NxSs),gridObj.y(NySs),cos(phiMaxSpat'), sin(phiMaxSpat'),...
            'color',[0 0 0],'AutoScaleFactor',0.5 );
        hold off
        TitlStr = sprintf('Min Max Plot t = %f N~%f', TimeRecVec(t), NumPart );
        title(TitlStr)
        
        subplot(3,1,2)
        pcolor(gridObj.x(NxSs),gridObj.y(NySs),rhoMin');
        set(gca,'CLim', [AbsMinRho AbsMaxRho],'YDir','rev');
        shading interp;
        colorbar;
        hold on
        quiver(gridObj.x(NxSs),gridObj.y(NySs),cos(phiMinSpat'), sin(phiMinSpat'),...
            'color',[0 0 0],'AutoScaleFactor',0.5 );
        hold off
        
        subplot(3,1,3)
        pcolor(gridObj.x(NxSs),gridObj.y(NySs), mean( Density_rec(NxSs,NySs,:,t),3 )' );
        set(gca,'CLim', [AbsMinRho AbsMaxRho],'YDir','rev');
        shading interp;
        colorbar;
        
        % Store Frame
        M_MinMaxSS(t) = getframe(gcf);
    end %end frame loop
    
end %MinMaxSSPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% keyboard

if AdjacentPointsPlot
    
    xPos1 = 17;
    yPos1 = 17;
    
    xPos2 = 18;
    yPos2 = 17;
    
    if MinMaxPlot == 0
        AbsMinRho = min( min( min( min( Density_rec ) ) ) );
        AbsMaxRho = max( max( max( max( Density_rec ) ) ) );
    end
    figure(3)
    set(gca,'NextPlot','replaceChildren',...
        'YLim',[ AbsMinRho AbsMaxRho ]);
    set(gcf,'renderer','zbuffer');
    
    for t = 1:nFrames
        
        plot( gridObj.phi, reshape( Density_rec(xPos1,yPos1,:,t), 1, systemObj.Nm ), ...
            gridObj.phi, reshape( Density_rec(xPos2,yPos2,:,t), 1, systemObj.Nm ) );
        
        TtlStr = sprintf('t = %f', TimeRecVec(t) );
        title(TtlStr);
        legStr1 = sprintf('xPos = %d yPos = %d', xPos1, yPos1);
        legStr2 = sprintf('xPos = %d yPos = %d', xPos2, yPos2);
        legend( legStr1 ,  legStr2);
        M_dist(t) = getframe(gcf);
        
    end % frame loop
    
end %AdjacentPointsPlot

end %function


