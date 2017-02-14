% [NegDivFluxEx_FT] = dRhoIntCalcVcFt(rho,muExFt,systemObj,diffObj)
%
% Calculate dRho given a chemical potentionial (FT) muExFt
%
function [NegDivFluxEx_FT] = ...
  dRhoIntCalcMu(rho,muExFt,systemObj,diffObj)
n3 = systemObj.n3;
% See if muExFt is only a position of position
[~,~,n3mu] = size(muExFt);
if n3mu > 1
  ik1 = diffObj.ik1rep3;
  ik2 = diffObj.ik2rep3;
  ik3 = diffObj.ik3rep3;
else
  ik1 = diffObj.ik1rep3(:,:,1);
  ik2 = diffObj.ik2rep3(:,:,1);
  ik3 = diffObj.ik3rep3(:,:,1);
end
%Takes its derivative in k-space, product in real, then back to k-space
% coordinate 1
dMuEx_dx1_FT   = ik1 .*  muExFt; % derivative of mu in k space
dMuEx_dx1   =  real(ifftn(ifftshift(dMuEx_dx1_FT))); % back to real
jx1 = - rho .* dMuEx_dx1;    %Flux in the x1 direction with isostropic diffusion
Jx1_FT = fftshift(fftn(jx1)); % FFT 
% coordinate 2
dMuEx_dx2_FT   = ik2 .*  muExFt;
dMuEx_dx2   =  real(ifftn(ifftshift(dMuEx_dx2_FT)));
jx2 = - rho .* dMuEx_dx2;    %Flux in the x2 direction with isostropic diffusion
Jx2_FT = fftshift(fftn(jx2));
% coordinate 3
if n3mu > 1
  dMuEx_dx3_FT = ik3 .*  muExFt;
  dMuEx_dx3 =  real(ifftn(ifftshift(dMuEx_dx3_FT)));
  jx3 = - rho .* dMuEx_dx3;  %Flux in the angular direction with isostropic diffusion
  Jx3_FT = fftshift(fftn(jx3));
end
% Do it all with indexing
if  n3 > 1 
  Ind    = [ 1:n3 ];
  Ind_m2 = [ n3-1, n3,  1:(n3-2) ]; %m-2 coupling
  Ind_p2 = [ 3:n3, 1, 2 ]; %m+2 coupling
else
  Ind = 1;
  Ind_m2 = 1;
  Ind_p2 = 1;
end
% Allocate and calulate
NegDivFluxEx_FT = zeros( systemObj.n1, systemObj.n2, n3mu );
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
if n3mu > 1
  NegDivFluxEx_FT = NegDivFluxEx_FT ...
    - ik3 .* diffObj.Mob_rot .* Jx3_FT;
end
