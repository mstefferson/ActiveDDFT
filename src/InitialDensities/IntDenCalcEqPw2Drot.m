% IntDen2DrotCalcPerturbEqPw.m 
% Description: creates the initial density in 2 spatial directions and one
% rotation DoF. The initial density of the form
% rho = rho_equilbrium + \sum A_k exp( i (k_x x + k_y y + k_phi phi)
% A_k is an input parameter


function [rho] = IntDenCalcEqPw2Drot(systemObj, rhoInit)

%Add in some slight deviation from the equilbrium density at specific modes.
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

% Build rho from equilibrium
% Initialize rho
rho = systemObj.c .* ...
    ones(systemObj.n1,systemObj.n2,systemObj.n3);

% Map distribution to a homogeneous system
for i = 1:systemObj.n3
    rho(:,:,i) = rho(:,:,i) .* rhoInit.feq(i);
end

end %end function
