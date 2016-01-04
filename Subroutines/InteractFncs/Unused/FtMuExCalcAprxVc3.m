% CURRENTLY UNUSED 4-Jan-2015
%
% Program: FtMuExCalcAprxVc
%
% Description: Heavily approximates the 3rd virial coefficient assuming that
% the density is constant
% 
% Called by: dRhoInterCalcVcID

function [MuEx_FT] = FtMuExCalcAprxVc3(rho,ParamObj)
%%%%%%%%%%%%%%%%%%%Hard rod interactions%%%%%%%%%%%%%%%%%%%%%%%%%%

% A drastic approximation of the excess chemical potential from the
% 3rd virial coefficient
MuEx_FT = - 1 / 2 .* (pi * 5 / 2 + 3) .* ParamObj.L_rod ^2 ...
    .* fftshift(fftn(rho .* rho));

end
