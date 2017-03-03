% [MuEx_FT] = muExCalcMfFt(rhoFt,vFt,systemObj,scaleFact)
%
% Description: Calculates the excess chemical potential using mean field
%
% Called by: dRhoMaster
% Excess chemical potential in position space is a convolution. In k-space, it is a
% product. Given by the function derivative of the excess free energy w.r.t.
% the density profile

function [MuEx_FT] = muExCalcMfFt(rhoFt,vFt,systemObj,scaleFact)
% Calc scale factor if need be
if nargin == 3
  scaleFact = (systemObj.l3 * systemObj.l1 * systemObj.l2) / (systemObj.n1 * systemObj.n2 * systemObj.n3);
end
MuEx_FT = scaleFact .* vFt .* rhoFt;
end
