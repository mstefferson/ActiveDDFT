% Function: MuExCalcVc2Ft.m
%
% Description: Calculates the excess chemical potential from Onsager's 2nd Virial Coefficient
% for a inhomogenous gas of hard rods.
%
% Called by: dRhoInterCalcVcID

function [MuEx_FT] = MuExCalcVc2Ft(rho_FT,Fm_FT,systemObj)


%Excess chemical potential in position space is a convolution. In k-space, it is a
%product. Given by the function derivative of the excess free energy w.r.t.
%the density profile
% keyboard
%Now includes the correct scale
MuEx_FT = -(systemObj.Lphi * systemObj.Lx * systemObj.Ly) / (systemObj.Nx * systemObj.Ny * systemObj.Nm) ...
    .* systemObj.Tmp .* Fm_FT .* rho_FT;

end
