% [MuEx_FT] = muExCalcVc2Ft(rho_FT,Fm_FT,systemObj,scaleFact)
%
% Description: Calculates the excess chemical potential 
% from direct correction function (mayer function expansion or mf)
% for a inhomogenous gas of hard rods.
%
% Called by: dRhoMaster
%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile

function [MuExFt] = muExCalcPDirCorrFt(rhoFt,c2Ft,systemObj,scaleFact)
%Now includes the correct scale
if nargin == 3
  scaleFact = (systemObj.l3 * systemObj.l1 * systemObj.l2) / (systemObj.n1 * systemObj.n2 * systemObj.n3);
end
MuExFt = - scaleFact .* c2Ft .* rhoFt;
end
