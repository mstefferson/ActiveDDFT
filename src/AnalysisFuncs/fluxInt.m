function [JxInt, JyInt, JphiInt] =  ...
  fluxInt( rho, rho_FT, D, diffObj, systemObj, particleObj )

% Mayer function stuff %
Fm_FT = fftshift(fftn( mayerFncHr(...
  systemObj.n1, systemObj.n2, systemObj.n3, ...
  systemObj.l1, systemObj.l2, particleObj.lMaj) ));

% Calc chemical potential
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

% Calculate isotropic "fluxes" then hit with diffusion matrix
jx = - rho .* dMuEx_dx;    %Flux in the x direction with isostropic diffusion
jy = - rho .* dMuEx_dy;    %Flux in the y direction with isostropic diffusion
jphi = - rho .* dMuEx_dphi;  %Flux in the angular direction with isostropic diffusio

% Now for flux
JxInt = ( D.xx .* jx + D.xy .* jy );
JyInt = ( D.xy .* jx + D.yy .* jy );
JphiInt =  ( D.mm .* jphi );

