% IntDen2DrotCalcSepPwModes.m
%
% Description: creates the initial density in 2 spatial directions and one
% rotation DoF. The initial density of the form
% rho = \sum A_k exp( i k_x x ) + \sum A_k exp( i k_y y ) ...
%       + \sum A_k exp( i k_m phi )
% A_k is an input parameter


function [rho] = IntDenCalcIsoSepPw2Drot(GridObj,ParamObj,RhoInit)

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

%Initialize rho
rho = ParamObj.Norm / (2 .* pi .* ParamObj.Lx .* ParamObj.Lx) .* ...
  ones(ParamObj.Nx,ParamObj.Ny,ParamObj.Nm);

% Perturb it
[rho] = SepPwDenPerturber2Drot(rho,ParamObj,GridObj,RhoInit);

% keyboard
end %end function
