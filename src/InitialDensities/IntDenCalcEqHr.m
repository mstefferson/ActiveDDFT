% IntDenCalcEqHr.m 
% Description: creates the initial density in 2 spatial directions and one
% rotation DoF. Angular distribution based on steady state Onsager theory

function [rho] = IntDenCalcEqHr(systemObj, rhoInit)
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
