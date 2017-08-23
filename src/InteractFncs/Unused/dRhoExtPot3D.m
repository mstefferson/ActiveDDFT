function [NegDivFluxPot_FT] = ...
  dRhoExtPot3D(rho, F, diffObj, systemObj)

%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile

jx = - rho .* F(:,:,:,1);    %Flux in the x direction with isostropic diffusion
jy = - rho .* F(:,:,:,2);    %Flux in the y direction with isostropic diffusion
jm = - rho .* F(:,:,:,3);  %Flux in the angular direction with isostropic diffusion

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

% Calc -grad j
NegDivFluxPot_FT = zeros( systemObj.n1, systemObj.n2, systemObj.n3 );
NegDivFluxPot_FT(:,:,Ind) = ...
  diffObj.j1f_reps .* Jx_FT(:,:,Ind) + ...
  diffObj.j1Mm2f_reps .* Jx_FT(:,:,Ind_m2) + ...
  diffObj.j1Mp2f_reps .* Jx_FT(:,:,Ind_p2) + ...
  diffObj.j2f_reps .* Jy_FT(:,:,Ind) + ...
  diffObj.j2Mm2f_reps .* Jy_FT(:,:,Ind_m2) + ...
  diffObj.j2Mp2f_reps .* Jy_FT(:,:,Ind_p2);

%Add the C(k) term last
NegDivFluxPot_FT = NegDivFluxPot_FT ...
  - diffObj.ik3rep3 .* diffObj.Mob_rot .* Jm_FT;

