% Function: MuExCalcVc2Ft.m
%
% Description: Calculates the excess chemical potential from Onsager's 2nd Virial Coefficient
% for a inhomogenous gas of hard rods.
%
% Called by: dRhoInterCalcVcID

function [MuEx_FT] = MuExCalcVc2Ft(rho_FT,Fm_FT,ParamObj)


%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile
% keyboard
%Now includes the correct scale
MuEx_FT = -(2 * pi * ParamObj.Lx * ParamObj.Ly) / (ParamObj.Nx * ParamObj.Ny * ParamObj.Nm) ...
    .* ParamObj.Tmp .* Fm_FT .* rho_FT;

end
