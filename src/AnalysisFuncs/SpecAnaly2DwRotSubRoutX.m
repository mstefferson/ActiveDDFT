function [ampl_record,fitParam_Re] = ...
    SpecAnaly2DwRotSubRoutX(NumModesMax,DensityFT_record,TimeRecVec, dt,kyholder,kmholder,...
    min_amp,kx,ky,km,D_rot, D_pos,Nx,Ny,Nm,bc)

DecayDisp   = 1;
AllKsVsTime = 1;

%Use squeeze to make a matrix (length(ky), length(j_record) ) of the
%amplitudes we want to look at
ampl_record = squeeze(DensityFT_record(:,kyholder,kmholder,:));

[k2plot,Nmodes] = Ks2TrackFinderSpec(NumModesMax,ampl_record,TimeRecVec,min_amp);

% keyboard
%Now plot all these amplitudes throughout the record


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DecayDisp
    
    % Get Dispersion Relations
%     lambda_c    = DispRelMakerContin(kx(k2plot),ky(kyholder),km(kmholder),D_pos,D_rot,bc);
%     lambda_dis  = DispRelMakerDiscrete(kx(k2plot),ky(kyholder),km(kmholder),D_pos,D_rot,bc,dt);
%     lambda_diff = DispRelMakerDiffus(kx(k2plot),ky(kyholder),km(kmholder),D_pos,D_rot);

lambda_diff = IsoDispRelMakerDiffus(kx(k2plot),ky(kyholder),km(kmholder),D_pos,D_rot);

lambda_c    = ...
    IsoDispRelMakerContin(kx(k2plot),ky(kyholder),km(kmholder),D_pos,D_rot,bc);
lambda_dis  = ...
    IsoDispRelMakerDiscrete(kx(k2plot),ky(kyholder),km(kmholder),D_pos,D_rot,bc,dt);


    
%     keyboard
    % Next figure is going to plot certain modes (same as
    % ampl_record(ky2plot)with respect to time. Along with this, we will plot a
    % exponential fit to each decaying mode.
    % Form of decay is: y = aexp(lambda*t)
    % where lambda = k^2 (1D)
    %define the fitting function and the parameters to save
    
    expfit = @(a,TimeRecVec) a(1)*exp(a(2)*TimeRecVec);  %Exponential fit function
    fitParam_Re = zeros(2,Nmodes);               %Parameter array
    options = optimset('Display','off'); %Don't display anything in fnc lsqcurvefit
    figure(1)
    subplot(2,1,1)
    
    for j=1:Nmodes
        % Make a vector of single mode amplitude throughout time.
        %Make them all positive so they all look similiar. Sign doesn't matter anyway
        %
        % Make a vector of single mode amplitude throughout time.
        %Make them all positive so they all look similiar. Sign doesn't matter anyway
        
        Re_yvec = abs(real(ampl_record(k2plot(j),:)));
        Re_a0 = [Re_yvec(1); -1 / TimeRecVec(end) * log( Re_yvec(1) / Re_yvec(end) )];
        %A guess for the fit parameters
        Coeff_fit_crnt_mode_Re = lsqcurvefit(expfit,Re_a0,TimeRecVec,Re_yvec,[],[],options);
        fitParam_Re(:,j) = Coeff_fit_crnt_mode_Re;  %Save the fit parameters
        
        if k2plot(j) ~= Nx/2 + 1
            plot(TimeRecVec,Re_yvec,'o-'); hold all;
            plot(TimeRecVec,Coeff_fit_crnt_mode_Re(1)*exp(Coeff_fit_crnt_mode_Re(2)*TimeRecVec),'-','LineWidth',2);
        end
        
    end %End Mode loop
    
    title('k_x mode amplitudes decay with fit, excluding k_x = 0');
    xlabel('t'); ylabel('Amplitude');
    hold off
    
    % Plot the kx = +1 modes seperately
    kxp1_holder = Nx/2+2;
    subplot(2,1,2)
    %     keyboard
    plot(TimeRecVec, abs(real( ampl_record(kxp1_holder ,:) )),'o', ...
        TimeRecVec,  abs(real( ampl_record(kxp1_holder ,1) )).* ...
        exp(lambda_c(k2plot == kxp1_holder) .* TimeRecVec) ,...
        TimeRecVec,  abs(real( ampl_record(kxp1_holder ,1) )).* ...
        exp(lambda_dis(k2plot == kxp1_holder) .* TimeRecVec)                 );
    legend('Measured','Predicted Contin LSA','Predicted Discrete LSA',...
        'Location','Best');
    titstr = sprintf('kx = %i mode vs time', kxp1_holder - (Nx /2 + 1) );
    title(titstr)
    xlabel('time')
    ylabel('amplitude')
    
    % keyboard
    % Plot the dispersion relation
    figure(2)
    plot(kx(k2plot),fitParam_Re(2,:),'o',kx(k2plot), lambda_c,'-', ...
        kx(k2plot), lambda_dis,'-',kx(k2plot), lambda_diff,'-') ;
    titstr = sprintf('k_x dispersion relation e^{lambda t}(k_y mode = %i, k_m mode = %i)',...
        kyholder - (Ny/2+1),kmholder - (Nm/2+1) );
    title(titstr)
    legend('Measured Disp Real','Predicted Continuous Linear Stab',...
        'Predicted Discrete Linear Stab','Just Diffustion Linear Stab',...
        'Location','Best');
    xlabel('k_x'); ylabel('Decay constant \lambda');
    %         keyboard
end % end DecayDisp

if AllKsVsTime
    
    figure(3)
    % Plot the mode amplitudes except for the k =
    NumTicks = 3; % Number of Tick marks on the y axis
    for j = 1:(Nmodes + 1 ) / 2
        % Make a vector of single mode amplitude throughout time.
        %Make them all positive so they all look similiar. Sign doesn't matter anyway
        
        if k2plot(j) == (Nx/2 + 1)
            
            subplot( (Nmodes+1) / 2, 2,[1,2] )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plot(j),:) ),...
                TimeRecVec,imag( ampl_record(k2plot(j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                k2plot(j) - ( Nx / 2 + 1 ), ...
                kyholder  - ( Ny / 2 + 1 ), ...
                kmholder  - ( Nm / 2 + 1 ) );
            ylabel(haxes(1),'Re \{F_k\}')
            ylabel(haxes(2),'Im \{F_k\}')
            xlabel(haxes(2),'time')
            set(haxes(1),'YTick', mean( real( ampl_record(k2plot(j),:) ) ) )
            set(haxes(2),'YTick',mean( imag( ampl_record(k2plot(j),:) ) ))
            title(titlestr)
            
        else
            
            %plot negative mode first
            subplot( (Nmodes+1) / 2, 2, (Nmodes + 3) -  2*j )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plot(j),:) ),...
                TimeRecVec,imag( ampl_record(k2plot(j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                k2plot(j) - ( Nx / 2 + 1 ), ...
                kyholder  - ( Ny / 2 + 1 ), ...
                kmholder  - ( Nm / 2 + 1 ) );
            ylabel(haxes(1),'Re \{F_k\}')
            ylabel(haxes(2),'Im \{F_k\}')
            %             xlabel(haxes(2),'time')
            WantedYTicks = linspace(...
                min( real( ampl_record(k2plot(j),:) ) ) , ...
                max( real( ampl_record(k2plot(j),:) ) ),  NumTicks);
            set(haxes(1),'YTick', WantedYTicks)
            
            WantedYTicks = linspace(...
                min( imag( ampl_record(k2plot(j),:) ) ) , ...
                max( imag( ampl_record(k2plot(j),:) ) ),  NumTicks);
            set(haxes(2),'YTick',WantedYTicks)
            title(titlestr)
            
            %Then plus
            subplot( (Nmodes+1) / 2, 2, (Nmodes + 2) -  2*j )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plot(Nmodes + 1 - j),:) ),...
                TimeRecVec,imag( ampl_record(k2plot(Nmodes + 1 - j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                k2plot(Nmodes + 1 - j) - ( Nx / 2 + 1 ), ...
                kyholder - ( Ny / 2 + 1 ), ...
                kmholder - ( Nm / 2 + 1 ) );
            ylabel(haxes(1),'Re \{F_k\}')
            ylabel(haxes(2),'Im \{F_k\}')
            %             xlabel(haxes(2),'time')
            
            
            WantedYTicks = linspace(...
                min( real( ampl_record(k2plot(Nmodes + 1 - j),:) ) ) , ...
                max( real( ampl_record(k2plot(Nmodes + 1 - j),:) ) ),  NumTicks);
            set(haxes(1),'YTick', WantedYTicks)
            
            WantedYTicks = linspace(...
                min( imag( ampl_record(k2plot(Nmodes + 1 - j),:) ) ) , ...
                max( imag( ampl_record(k2plot(Nmodes + 1 - j),:) ) ),  NumTicks);
            set(haxes(2),'YTick',WantedYTicks)
            
            title(titlestr)
        end %end if k =0 statement
    end %end mode loop
    
end % end AllKsVsTime

% keyboard
end %end function