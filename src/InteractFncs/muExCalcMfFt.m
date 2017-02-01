% Function: MuExCalcVc2Ft.m
%
% Description: Calculates the excess chemical potential from Onsager's 2nd Virial Coefficient
% for a inhomogenous gas of hard rods.
%
% Called by: dRhoInterCalcVcID
%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile

function [MuEx_FT] = muExCalcMfFt(rhoFt,vFt,systemObj,scaleFact)
% Calc scale factor if need be
if nargin == 3
  scaleFact = (systemObj.Lphi * systemObj.Lx * systemObj.Ly) / (systemObj.Nx * systemObj.Ny * systemObj.Nm);
end
MuEx_FT = scaleFact .* vFt .* rhoFt;
end
