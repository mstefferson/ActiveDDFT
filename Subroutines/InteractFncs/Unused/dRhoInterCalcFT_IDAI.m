%% CURRENTLY NOT USED 4-Jan-16

% Function: dRhoInterCalcFT_AI.m

% Description: Uses the 2nd virial coefficient to calcuate the change in
% density due to hard rod interactions in k-space assuming spatial homogenity.
% This should only be used if there are no spatial variation. Otherwise, the
% assumptions could lead to issues. Uses traditional Onsager excluded volume:
% kernal is sin( d_phi ) instead of Mayer function.
% 
% ISOTROPIC diffusion

function [NegDivFluxExcess_FT] = dRhoInterCalcFT_IDAI(rho_FT,ParamObj,GridObj,DiffMobObj)
%%%%%%%%%%%%%%%%%%%Hard rod interactions%%%%%%%%%%%%%%%%%%%%%%%%%%

%Excess chemical potential in orientational space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile

% Calculate rho
rho = real( ifftn(ifftshift( rho_FT ) ) );

% Just integrating over angle. Already integrated over space to get kernal.
% Just Fourier transform w.r.t. to angle
% keyboard
Kern_FT  = fftshift( fft( abs( sin( GridObj.phi3D) ),[], 3) );

%Now includes the correct scale
MuEx_FT = (2 * pi ) / ParamObj.Nm .* ParamObj.Tmp .* ParamObj.L_rod ^ 2 ...
    .* Kern_FT .* rho_FT;

%     MuEx    = real(ifftn(ifftshift(MuEx_FT)));

%Takes its derivative in k-space
dMuEx_dx_FT   =     sqrt(-1) .* GridObj.kx3D .*  MuEx_FT;
dMuEx_dy_FT   =     sqrt(-1) .* GridObj.ky3D .*  MuEx_FT;
dMuEx_dphi_FT =     sqrt(-1) .* GridObj.km3D .*  MuEx_FT;

%Excess chemical potential derivative in real space
%Mayer function derivative in real-space
dMuEx_dx   =  real(ifftn(ifftshift(dMuEx_dx_FT)));
dMuEx_dy   =  real(ifftn(ifftshift(dMuEx_dy_FT)));
dMuEx_dphi =  real(ifftn(ifftshift(dMuEx_dphi_FT)));

%Do the hard disk interaction portion of the PDE in real space then FT it
% Isolate the seperate parts and call them some arbitrary function. We
% will Fourier transform these functions to solve this in Fourier space
%
% Take the divergence of the product of functions. Call these products
% random variables

jx = - DiffMobObj.Mob_pos   .* rho .* dMuEx_dx;    %Flux in the x direction with isostropic diffusion
jy = - DiffMobObj.Mob_pos   .* rho .* dMuEx_dy;    %Flux in the y direction with isostropic diffusion
jm = - DiffMobObj.Mob_rot   .* rho .* dMuEx_dphi;   %Flux in the angular direction with isostropic diffusion

%Fourier transform these
Jx_FT = fftshift(fftn(jx));
Jy_FT = fftshift(fftn(jy));
Jm_FT = fftshift(fftn(jm));

% Calculate the - divergence of the interaction flux
NegDivFluxExcess_FT = - sqrt(-1) .* ( GridObj.kx3D .* Jx_FT + ...
    GridObj.ky3D .* Jy_FT + GridObj.km3D .* Jm_FT );

