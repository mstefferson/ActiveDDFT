function [rho] = IntDenCalcGauss2Drot(gridObj,systemObj,rhoInit)

%Add in some slight deviation from a uniform density at specific modes.
% Here, our perturbation is a gaussian.
% Now for modes: 0 means uniform, 1 means Gaussian
% Density is normalized so that
%
% # of particles           = int( rho(x,y,phi) dx dy dphi )
% c(x,y) (concentration)   = int( rho(x,y,phi) dphi )
% f(phi) (AngDistribution) = int( rho(x,y,phi) dx dy ) ./ # Particles
% 1                        = int( f(phi) dphi)

var = .1;

rho = ones(systemObj.Nx,systemObj.Ny,systemObj.Nm); % k = 0

% Pick a random angle to center Gaussian on
AngMax = 2*pi * rand();

%find closest angle
ClosestMaxAngInd = find( gridObj.phi  - AngMax ==  min(abs( gridObj.phi - AngMax) ) ) ;
GaussVecTemp = [ exp( - ( gridObj.phi(1:systemObj.Nm/2) ) .^2 / var^2 ) ...
  exp( - ( gridObj.phi(systemObj.Nm/2+1:end) - 2*pi ) .^2 / var^2 )];
distTemp = circshift(GaussVecTemp',ClosestMaxAngInd)';

for i = 1:systemObj.Nm
  rho(:,:,i) = distTemp(i);
end

% keyboard
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(gridObj.y,trapz_periodic(gridObj.x,trapz_periodic(gridObj.phi,rho,3),2),1);
rho = rho .* systemObj.numPart ./ CurrentNorm;


end %end function
