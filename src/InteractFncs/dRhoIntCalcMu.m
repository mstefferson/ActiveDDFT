% [NegDivFluxEx_FT] = dRhoIntCalcVcFt(rho,muExFt,systemObj,diffObj)
%
% Calculate dRho given a chemical potentionial (FT) muExFt
%
function [negDivFluxEx_FT, jx3] = ...
  dRhoIntCalcMu(rho, dMu, systemObj, diffObj, interObj)
% Allocate and calulate
n3 = systemObj.n3;
% bug dMu does not have a size
negDivFluxEx_FT = zeros( systemObj.n1, systemObj.n2, n3 );
% Do it all with indexing
if  n3 > 1
  jInd    = [ 1:n3 ];
  jInd_m2 = [ n3-1, n3,  1:(n3-2) ]; %m-2 coupling
  jInd_p2 = [ 3:n3, 1, 2 ]; %m+2 coupling
else
  jInd = 1;
  jInd_m2 = 1;
  jInd_p2 = 1;
end
% coordinate 1
if interObj.dv1Flag
  jx1 = - rho .* dMu.dx1;    %Flux in the x1 direction with isostropic diffusion
  jx1_FT = fftshift(fftn(jx1)); % FFT
  if diffObj.Ani == 0
    negDivFluxEx_FT(:,:,jInd) = ...
      diffObj.j1f_reps .* jx1_FT(:,:,jInd);
  else
    negDivFluxEx_FT(:,:,jInd) = ...
      diffObj.j1f_reps .* jx1_FT(:,:,jInd) + ...
      diffObj.j1Mm2f_reps .* jx1_FT(:,:,jInd_m2) + ...
      diffObj.j1Mp2f_reps .* jx1_FT(:,:,jInd_p2);
  end
end
% coordinate 2
if interObj.dv2Flag
  jx2 = - rho .* dMu.dx2;    %Flux in the x2 direction with isostropic diffusion
  jx2_FT = fftshift(fftn(jx2));
  if diffObj.Ani == 0
    negDivFluxEx_FT(:,:,jInd) = negDivFluxEx_FT(:,:,jInd) + ...
      diffObj.j2f_reps .* jx2_FT(:,:,jInd);
  else
    negDivFluxEx_FT(:,:,jInd) = negDivFluxEx_FT(:,:,jInd) +...
      diffObj.j2f_reps .* jx2_FT(:,:,jInd) + ...
      diffObj.j2Mm2f_reps .* jx2_FT(:,:,jInd_m2) + ...
      diffObj.j2Mp2f_reps .* jx2_FT(:,:,jInd_p2);
  end
end
% coordinate 3
if interObj.dv3Flag
  jx3 = - rho .* dMu.dx3;  %Flux in the angular direction with isostropic diffusion
  jx3_FT = fftshift(fftn(jx3));
  negDivFluxEx_FT = negDivFluxEx_FT +...
    diffObj.j3f_reps .* jx3_FT;
end
