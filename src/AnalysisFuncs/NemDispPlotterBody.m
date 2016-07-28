function NemDispPlotterBody(GridObj,ParamObj,k2plotInd,ampl_record,Nmodes,TimeRecVec, kxholder,kyholder,...
    kx,ky,km,D_rot, D_pos,Nx,Ny,Nm,bc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     keyboard
[lambda_c] = ...
    NemDispRelMakerContin(GridObj,ParamObj,ampl_record,kx(kxholder),ky(kyholder),km(k2plotInd),...
    km,Nm,D_pos,D_rot,bc);
[lambda_diff] = IsoDispRelMakerDiffus(kx(kxholder),ky(kyholder),km(k2plotInd),D_pos,D_rot);

% Next figure is going to plot certain modes (same as
% ampl_record(ky2plot)with respect to time. Along with this, we will plot a
% exponential fit to each decaying mode.
% Form of decay is: y = aexp(lambda*t)
% where lambda = k^2 (1D)
%define the fitting function and the parameters to save
%Exponential fit function; epsilon = a(1) + i a(2); lambda = a(3) + i a(4);
Re_expfit = @(a,TimeRecVec) ...
    exp(a(3).*TimeRecVec) .* ( a(1).* cos( a(4) .* TimeRecVec )  ...
    - a(2).*sin( a(4) .* TimeRecVec ) ) ;
Im_expfit = @(a,TimeRecVec) ...
    exp(a(3).*TimeRecVec) .* ( a(2).*cos( a(4) .* TimeRecVec )  ...
    + a(1).*sin( a(4) .* TimeRecVec ) ) ;

fitParam_Re = zeros(4,Nmodes);               %Parameter array
options = optimset('Display','off'); %Don't display anything in fnc lsqcurvefit
figure(7)

% Need to track deviation from equilbrium, just, equilibrium
[rho_eq] = DenEq2Drot(GridObj,ParamObj);
kx_holder = 17;
ky_holder = 17;

FT_rho_eq = fftshift( fftn( rho_eq ) );
FT_rho_eqVec = reshape(  FT_rho_eq(kx_holder,ky_holder,:), 1, 32 );

% epsilon = FT_orig - FT_rho_eqVec;

subplot(1,1,1)

for j=1:Nmodes
    % Make a vector of single mode amplitude throughout time.
    %Make them all positive so they all look similiar. Sign doesn't matter anyway
    %
    % Make a vector of single mode amplitude throughout time.
    %Make them all positive so they all look similiar. Sign doesn't matter anyway
    
    yvec_t = ampl_record(k2plotInd(j),:);
    
    %epsilon(t)
    eps_vec_t = yvec_t - FT_rho_eqVec( k2plotInd(j) );
%     yvec_t(1) - FT_rho_eqVec( k2plotInd(j) );
    % Make a guess for the fit
    
%     keyboard
%     Re_a0 = [ real( yvec(1) ) ; ...
%         imag( yvec(1) )  ; ...
%         real( lambda_c(j) ) ; ...
%         imag( lambda_c(j) )   ];
%     

 Re_a0 = [ real( eps_vec_t(1) ) ; ...
        imag( eps_vec_t(1) )  ; ...
        real( lambda_c(j) ) ; ...
        imag( lambda_c(j) )   ];
   
    %A guess for the fit parameters
    Coeff_fit_crnt_mode_Re = lsqcurvefit(Re_expfit,Re_a0,TimeRecVec,real(eps_vec_t),[],[],options);
    fitParam_Re(:,j) = Coeff_fit_crnt_mode_Re;  %Save the fit parameters
    
%     keyboard
    if k2plotInd(j) ~= Nm/2 + 1
        plot(TimeRecVec,real(eps_vec_t),'o-',TimeRecVec,real(yvec_t),'x-');
        hold all;
        plot(TimeRecVec, ...
            exp(Coeff_fit_crnt_mode_Re(3) .* TimeRecVec) .* ...
          ( Coeff_fit_crnt_mode_Re(1) .* cos( Coeff_fit_crnt_mode_Re(4) .* TimeRecVec ) ...
          - Coeff_fit_crnt_mode_Re(2) .* sin( Coeff_fit_crnt_mode_Re(4) .* TimeRecVec ) )...
                ,'-','LineWidth',2 );
        titstr = sprintf('km = %i',km( k2plotInd(j) ) );
        title(titstr);
        hold off
%         keyboard
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
% kmp2_holder = Nm/2+3;
% subplot(2,1,2)
% %     keyboard
% plot(TimeRecVec, abs(real( ampl_record(kmp2_holder ,:) )),'o', ...
%     TimeRecVec,  abs(real( ampl_record(kmp2_holder ,1) )).* ...
%     exp( lambda_c(k2plotInd == kmp2_holder) .* TimeRecVec ) ,...
%     TimeRecVec,  abs(real( ampl_record(kmp2_holder ,1) )).* ...
%     exp(lambda_dis(k2plotInd == kmp2_holder) .* TimeRecVec)                 );
% legend('Measured','Predicted Contin LSA','Predicted Discrete LSA',...
%     'Location','Best');
% titstr = sprintf('km = %i mode vs time', kx( kmp2_holder ) - (Nm /2 + 1) );
% title(titstr)
% title('km = 2 mode vs time')
% xlabel('time')
% ylabel('amplitude')

% keyboard
% Plot the dispersion relation
figure(8)
subplot(2,1,1)
plot(km(k2plotInd),fitParam_Re(3,:),'o',km(k2plotInd), real(lambda_c),'-', ...
    km(k2plotInd), lambda_diff,'-') ;
titstr = sprintf('Real lambda k_m dispersion relation e^{lambda t}(k_x mode = %i, k_y mode = %i)',...
    kxholder - (Nx/2+1),kyholder - (Ny/2+1) );
title(titstr)
legend('Real: Measured Disp','Real: Predicted Con. Linear Stab',...
    'Diffusion','Location','Best');
xlabel('k_m'); ylabel('Decay constant \lambda');

subplot(2,1,2)
plot(km(k2plotInd),fitParam_Re(4,:),'o',km(k2plotInd), imag(lambda_c),'-');
titstr = sprintf('Imag lambda k_m dispersion relation e^{lambda t}(k_x mode = %i, k_y mode = %i)',...
    kxholder - (Nx/2+1),kyholder - (Ny/2+1) );
title(titstr)
legend('Imag: Measured Disp','Imag: Predicted Con. Linear Stab',...
      'Location','Best');
xlabel('k_m'); ylabel('Decay constant \lambda');

%         keyboard

end % end function

