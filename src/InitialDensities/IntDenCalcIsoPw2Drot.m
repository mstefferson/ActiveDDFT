% IntDen2DrotCalcPwModes.m
%
% Description: creates the initial density in 2 spatial directions and one
% rotation DoF. The initial density of the form
% rho = \sum A_k exp( i (k_x x + k_y y + k_phi phi)
% A_k is an input parameter

function [rho] = IntDenCalcIsoPw2Drot(systemObj)

%Add in some slight deviation from a uniform density at specific modes.
% The number of modes counts the modes above and below k=0. But given the
% symmetry, these modes are the same if you add a perturbation like cos(kx)
% But really, we are adding 2*NumModes to the system
%
% Density is normalized so that
%
% # of particles           = int( rho(x,y,phi) dx dy dphi )
% c(x,y) (concentration)   = int( rho(x,y,phi) dphi )
% f(phi) (AngDistribution) = int( rho(x,y,phi) dx dy ) ./ # Particles
% 1                        = int( f(phi) dphi)

% Initial rho
rho = systemObj.numPart / (systemObj.l3 .* systemObj.l1 .* systemObj.l1) .* ...
          ones(systemObj.n1,systemObj.n2,systemObj.n3); % k = 0
%Recall that the k=0 mode is located at points N/2 + 1


end %end function
