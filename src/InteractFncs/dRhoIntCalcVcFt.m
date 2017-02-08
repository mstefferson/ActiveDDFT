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
  dRhoIntCalcVcFt(rho,rho_FT,Fm_FT,systemObj,diffObj,muScale)
n3 = systemObj.n3;
%%%%%%%%%%%%%%%%%%%Hard rod interactions%%%%%%%%%%%%%%%%%%%%%%%%%%
%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile
%Now includes the correct scale
MuEx_FT = muExCalcVc2Ft(rho_FT, Fm_FT,systemObj,muScale);
%Takes its derivative in k-space
dMuEx_dx1_FT   = diffObj.ik1rep3 .*  MuEx_FT;
dMuEx_dx2_FT   = diffObj.ik2rep3 .*  MuEx_FT;
dMuEx_dx3_FT = diffObj.ik3rep3 .*  MuEx_FT;
%Excess chemical potential derivative in real space
%Mayer function derivative in real-space
dMuEx_dx1   =  real(ifftn(ifftshift(dMuEx_dx1_FT)));
dMuEx_dx2   =  real(ifftn(ifftshift(dMuEx_dx2_FT)));
dMuEx_dx3 =  real(ifftn(ifftshift(dMuEx_dx3_FT)));
%Do the hard disk interaction portion of the PDE in real space then FT.
% Isolate the seperate parts and call them some arbitrary function. We
% will Fourier transform these functions to solve this in Fourier space
%
% Take the divergence of the product of functions. Call these products
% random variables
jx1 = - rho .* dMuEx_dx1;    %Flux in the x direction with isostropic diffusion
jx2 = - rho .* dMuEx_dx2;    %Flux in the y direction with isostropic diffusion
jx3 = - rho .* dMuEx_dx3;  %Flux in the angular direction with isostropic diffusion
%Fourier transform these
Jx1_FT = fftshift(fftn(jx1));
Jx2_FT = fftshift(fftn(jx2));
Jx3_FT = fftshift(fftn(jx3));
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
if diffObj.Ani == 0
  NegDivFluxEx_FT(:,:,Ind) = ...
    diffObj.j1f_reps .* Jx1_FT(:,:,Ind) + ...
    diffObj.j2f_reps .* Jx2_FT(:,:,Ind);
else
  NegDivFluxEx_FT(:,:,Ind) = ...
    diffObj.j1f_reps .* Jx1_FT(:,:,Ind) + ...
    diffObj.j1Mm2f_reps .* Jx1_FT(:,:,Ind_m2) + ...
    diffObj.j1Mp2f_reps .* Jx1_FT(:,:,Ind_p2) + ...
    diffObj.j2f_reps .* Jx2_FT(:,:,Ind) + ...
    diffObj.j2Mm2f_reps .* Jx2_FT(:,:,Ind_m2) + ...
    diffObj.j2Mp2f_reps .* Jx2_FT(:,:,Ind_p2);
end
%Add the C(k) term last
NegDivFluxEx_FT = NegDivFluxEx_FT ...
  - diffObj.ik3rep3 .* diffObj.Mob_rot .* Jx3_FT;
