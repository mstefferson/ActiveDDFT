% Function: MuExDisksSPT.m
%
% Description: Calculates the excess chemical potential for hard disks from scaled particle theory


function [MuEx_FT] = MuExDisksSPT(rho)

% Scaled particle theory
MuEx = - log( 1 - rho ) + ( 3 * rho - 2 * rho .^ 2 ) ./ ( 1 - rho ) ;
% Fourier
MuEx_FT = fftshift( fftn ( MuEx ) );

end
