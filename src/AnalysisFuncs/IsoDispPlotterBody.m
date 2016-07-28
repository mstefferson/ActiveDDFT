function IsoDispPlotterBody(k2plotInd,ampl_record,Nmodes,TimeRecVec, dt,kxholder,kyholder,...
    kx,ky,km,D_rot, D_pos,Nx,Ny,Nm,bc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get Dispersion Relations
lambda_diff = IsoDispRelMakerDiffus(kx(kxholder),ky(kyholder),km(k2plotInd),D_pos,D_rot);

lambda_c    = ...
    IsoDispRelMakerContin(kx(kxholder),ky(kyholder),km(k2plotInd),D_pos,D_rot,bc);
lambda_dis  = ...
    IsoDispRelMakerDiscrete(kx(kxholder),ky(kyholder),km(k2plotInd),D_pos,D_rot,bc,dt);


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
figure(7)
subplot(2,1,1)

for j=1:Nmodes
    % Make a vector of single mode amplitude throughout time.
    %Make them all positive so they all look similiar. Sign doesn't matter anyway
    %
    % Make a vector of single mode amplitude throughout time.
    %Make them all positive so they all look similiar. Sign doesn't matter anyway
    
    Re_yvec = abs(real(ampl_record(k2plotInd(j),:)));
    % Make a for the decaying/growing exponential
    %         if Re_yvec(1) > Re_yvec(end) %Decaying
    Re_a0 = [Re_yvec(1); -1 / TimeRecVec(end) * log( Re_yvec(1) / Re_yvec(end) )];
    %         else %Growing
    
    %         end
    
    %         keyboard
    %A guess for the fit parameters
    Coeff_fit_crnt_mode_Re = lsqcurvefit(expfit,Re_a0,TimeRecVec,Re_yvec,[],[],options);
    fitParam_Re(:,j) = Coeff_fit_crnt_mode_Re;  %Save the fit parameters
    
    if k2plotInd(j) ~= Nm/2 + 1
        plot(TimeRecVec,Re_yvec,'o-');
        hold all;
        plot(TimeRecVec,Coeff_fit_crnt_mode_Re(1)*exp(Coeff_fit_crnt_mode_Re(2)*TimeRecVec),'-','LineWidth',2);
        hold off
    end
    %         k2plot(j) - (Nm/2 + 1)
    %         Coeff_fit_crnt_mode_Re(2)
    %         Re_a0
    %         keyboard
end %End Mode loop

title('k_m mode amplitudes decay with fit, excluding k_m = 0');
xlabel('t'); ylabel('Amplitude');
hold off

%     keyboard
% Plot an individual mode
kmp2_holder = Nm/2+3;
subplot(2,1,2)
%     keyboard
plot(TimeRecVec, abs(real( ampl_record(kmp2_holder ,:) )),'o', ...
    TimeRecVec,  abs(real( ampl_record(kmp2_holder ,1) )).* ...
    exp( lambda_c(k2plotInd == kmp2_holder) .* TimeRecVec ) ,...
    TimeRecVec,  abs(real( ampl_record(kmp2_holder ,1) )).* ...
    exp(lambda_dis(k2plotInd == kmp2_holder) .* TimeRecVec)                 );
legend('Measured','Predicted Contin LSA','Predicted Discrete LSA',...
    'Location','Best');
titstr = sprintf('km = %i mode vs time', kx( kmp2_holder ) - (Nm /2 + 1) );
title(titstr)
title('km = 2 mode vs time')
xlabel('time')
ylabel('amplitude')

% keyboard
% Plot the dispersion relation
figure(8)
plot(km(k2plotInd),fitParam_Re(2,:),'o',km(k2plotInd), lambda_c,'-', ...
    km(k2plotInd), lambda_dis,'-',km(k2plotInd), lambda_diff,'-') ;
titstr = sprintf('k_m dispersion relation e^{lambda t}(k_x mode = %i, k_y mode = %i)',...
    kxholder - (Nx/2+1),kyholder - (Ny/2+1) );
title(titstr)
legend('Measured Disp Real','Predicted Continuous Linear Stab',...
    'Predicted Discrete Linear Stab','Just Diffustion Linear Stab',...
    'Location','Best');
xlabel('k_m'); ylabel('Decay constant \lambda');
%         keyboard

end % end function

