% Function: dRhoIntCalcVcFt.

% Description: Uses the 2nd virial coefficient to calcuate the change i
% density due to hard rod interactions in k-space.
%
% ANISOTROPIC Diffusion
%
% Called by: Anisotropic main
%
% Calls: MuExCalcVc2Ft

function [NegDivFluxEx_FT] = ...
  dRhoIntCalcVcFt(rho,rho_FT,Fm_FT,systemObj,diffObj)
n3 = systemObj.n3;
%%%%%%%%%%%%%%%%%%%Hard rod interactions%%%%%%%%%%%%%%%%%%%%%%%%%%

%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile

%Now includes the correct scale
MuEx_FT = muExCalcVc2Ft(rho_FT, Fm_FT,systemObj);

%Takes its derivative in k-space
dMuEx_dx_FT   = diffObj.ik1rep3 .*  MuEx_FT;
dMuEx_dy_FT   = diffObj.ik2rep3 .*  MuEx_FT;
dMuEx_dphi_FT = diffObj.ik3rep3 .*  MuEx_FT;

%Excess chemical potential derivative in real space
%Mayer function derivative in real-space
dMuEx_dx   =  real(ifftn(ifftshift(dMuEx_dx_FT)));
dMuEx_dy   =  real(ifftn(ifftshift(dMuEx_dy_FT)));
dMuEx_dphi =  real(ifftn(ifftshift(dMuEx_dphi_FT)));

%Do the hard disk interaction portion of the PDE in real space then FT.
% Isolate the seperate parts and call them some arbitrary function. We
% will Fourier transform these functions to solve this in Fourier space
%
% Take the divergence of the product of functions. Call these products
% random variables

jx = - rho .* dMuEx_dx;    %Flux in the x direction with isostropic diffusion
jy = - rho .* dMuEx_dy;    %Flux in the y direction with isostropic diffusion
jm = - rho .* dMuEx_dphi;  %Flux in the angular direction with isostropic diffusion

%Fourier transform these
Jx_FT = fftshift(fftn(jx));
Jy_FT = fftshift(fftn(jy));
Jm_FT = fftshift(fftn(jm));

%Calculate the -diverance of the flux in Fourier space. ;
% Do it all with indexing
if n3 > 1
  Ind    = [ 1:n3 ];
  Ind_m2 = [ n3-1, n3,  1:(n3-2) ]; %m-2 coupling
  Ind_p2 = [ 3:n3, 1, 2 ]; %m+2 coupling
else
  Ind = 1;
  Ind_m2 = 1;
  Ind_p2 = 1;
end
NegDivFluxEx_FT = zeros( systemObj.n1, systemObj.n2, n3 );
NegDivFluxEx_FT(:,:,Ind) = ...
  diffObj.jxf_reps .* Jx_FT(:,:,Ind) + ...
  diffObj.jxMm2f_reps .* Jx_FT(:,:,Ind_m2) + ...
  diffObj.jxMp2f_reps .* Jx_FT(:,:,Ind_p2) + ...
  diffObj.jyf_reps .* Jy_FT(:,:,Ind) + ...
  diffObj.jyMm2f_reps .* Jy_FT(:,:,Ind_m2) + ...
  diffObj.jyMp2f_reps .* Jy_FT(:,:,Ind_p2);

%Add the C(k) term last
NegDivFluxEx_FT = NegDivFluxEx_FT ...
  - diffObj.ik3rep3 .* diffObj.Mob_rot .* Jm_FT;

