% [NegDivFluxEx_FT] = dRhoIntCalcVcFt(rho,muExFt,systemObj,diffObj)
%
% Calculate dRho given iota where j = mob * iota.
% iota: the flux before getting hit by the mobility tensor, 

function [negDivFluxEx_FT] = ...
  dRhoFromFlux(rho, fluxFlag, diffObj, iota1, iota2, iota3)
% Allocate and calulate
n3 = systemObj.n3;
% bug dMu does not have a size
negDivFluxEx_FT = zeros( systemObj.n1, systemObj.n2, n3 );
% calc j1
if flag.calcj1
  % handle mixing if there is any
  if diffObj.Ani == 0
    j1 = diffObj.mob11 .* iota1 + diffObj.mob12 .* iota2;
  else
    j1 = diffObj.mob11 .* iota1;
  end
  j1Ft = fftshift( fftn( j1 ) );
  % calculate 
  negDivFluxEx_FT = - ( diffObj.ik1 .* j1Ft );
end
% calc j1
if flag.calcj2
  % handle mixing if there is any
  if diffObj.Ani == 0
    j2 = diffObj.mob12 .* iota1 + diffObj.mob22 .* iota2;
  else
    j2 = diffObj.mob22 .* iota2;
  end
  j2Ft = fftshift( fftn( j2 ) );
  % calculate 
  negDivFluxEx_FT = negDivFluxEx_FT - ( diffObj.ik2 .* j2Ft );
end
% calc j3
if flag.calcj3
  j3 = diffObj.mob33 .* iota3;
  j3Ft = fftshift( fftn( j3 ) );
  % calculate 
  negDivFluxEx_FT = negDivFluxEx_FT - ( diffObj.ik3 .* j3Ft );
end
