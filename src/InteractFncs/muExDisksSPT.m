% Function: MuExDisksSPT.m
%
% Description: Calculates the excess chemical potential for hard disks from scaled particle theory


function [muExFt] = muExDisksSPT(nu)
% Scaled particle theory
muEx = - log( 1 - nu ) + ( 3 * nu - 2 * nu .^ 2 ) ./ ( 1 - nu ) .^ 2 ;
% Fourier
muExFt = fftshift( fftn ( muEx ) );

end
