% [NegDivFluxEx_FT] = dRhoIntCalcVcFt(rho,muExFt,systemObj,diffObj)
%
% Calculate j given a chemical potential gradient
%
function [iota1, iota2, iota3] = ...
  fluxMu(rho, dMu,interObj)
% coordinate 1
if interObj.dv1Flag
  iota1 = - rho .* dMu.dx1; %Flux in the x1 direction before mobility
else
  iota1 = 0;
end
% coordinate 2
if interObj.dv2Flag
  iota2 = - rho .* dMu.dx2; %Flux in the x1 direction before mobility
else
  iota2 = 0;
end
% coordinate 3
if interObj.dv3Flag
  iota3 = - rho .* dMu.dx3; %Flux in the x1 direction before mobility
else
  iota3 = 0;
end
