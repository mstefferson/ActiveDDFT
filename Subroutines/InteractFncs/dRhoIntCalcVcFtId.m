% Function: dRhoInterCalcVcId.m

% Description: Wrapper function to calcuate the interaction contribution
% to the PDE in k-space. This calls functions that calculate the excess 
% chemical potential calculated by the virial expansion. It can call
% multiple expansions, but it should not be trusted above 2nd order. From
% excess chemical potential, it calculates the flux, then take the
% divergence of it, assuming isotropic diffusion.
% 
% ISOTROPIC Diffusion
% 
% Called by:  HR2DrotDenEvolverFTBodyIDCube.m - main body of isotropic
% diffusion

% Calls: 
% MuExCalcVc2Id(rho_FT, Fm_FT, ParamObj)- Calculates excess Chem pot. Using
% 2nd virial approx.

function [NegDivFluxExcess_FT] = ...
       dRhoIntCalcVcFtId(rho,rho_FT,Fm_FT,ParamObj,GridObj,DiffMobObj)

%%%%%%%%%%%%%%%%%%%Hard rod %%%%%%%%%%%%%%%%

%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile
% keyboard
%Now includes the correct scale

[MuEx_FT] = MuExCalcVc2Ft(rho_FT,Fm_FT,ParamObj);
% [MuEx3_FT] = FtMuExCalcAprxVc3(rho,ParamObj);
% [MuEx3_FT] = 0;
% keyboard

% MuEx_FT    = MuEx2_FT+MuEx3_FT;
%Takes its derivative in k-space
dMuEx_dx_FT   = DiffMobObj.ikx3 .*  MuEx_FT;
dMuEx_dy_FT   = DiffMobObj.iky3 .*  MuEx_FT;
dMuEx_dphi_FT = DiffMobObj.ikm3 .*  MuEx_FT;

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

jx = - DiffMobObj.Mob_pos .* rho .* dMuEx_dx;    %Flux in the x direction with isostropic diffusion
jy = - DiffMobObj.Mob_pos .* rho .* dMuEx_dy;    %Flux in the y direction with isostropic diffusion
jm = - DiffMobObj.Mob_rot .* rho .* dMuEx_dphi;  %Flux in the angular direction with isostropic diffusion

%Fourier transform these
Jx_FT = fftshift(fftn(jx));
Jy_FT = fftshift(fftn(jy));
Jm_FT = fftshift(fftn(jm));

% Calculate the - divergence of the interaction flux
NegDivFluxExcess_FT = - ( DiffMobObj.ikx3 .* Jx_FT + ...
    DiffMobObj.iky3 .* Jy_FT + DiffMobObj.ikm3 .* Jm_FT );

