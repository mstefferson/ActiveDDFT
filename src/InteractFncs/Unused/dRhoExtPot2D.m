function [NegDivFluxPot_FT] = ...
  dRhoExtPot3D(rho, F, diffObj, systemObj)

%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile

jx = - diffObj.Mob_pos .* rho .* F(:,:,1);    %Flux in the x direction with isostropic diffusion
jy = - diffObj.Mob_pos .* rho .* F(:,:,2);    %Flux in the y direction with isostropic diffusion

%Fourier transform these
Jx_FT = fftshift(fftn(jx));
Jy_FT = fftshift(fftn(jy));

%Calculate the -diverance of the flux in Fourier space. ;
NegDivFluxPot_FT = ...
  diffObj.ikx2 .* Jx_FT + ...
  diffObj.iky2 .* Jy_FT;
