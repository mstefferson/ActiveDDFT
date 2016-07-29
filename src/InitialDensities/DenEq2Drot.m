% IntDen2DrotCalcPerturbEqPw.m 
% Description: creates the initial density in 2 spatial directions and one
% rotation DoF. The initial density of the form
% rho = rho_equilbrium + \sum A_k exp( i (k_x x + k_y y + k_phi phi)
% A_k is an input parameter


function [rho_eq] = DenEq2Drot(gridObj,ParamObj)

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

% Add a path
addpath('C:\Users\MWS\Documents\MATLAB\Research\BG\INtransEq')
% Distribution stuff
Nc    = 20;            % Number of Coefficients

[Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,gridObj.phi,systemObj.bc); % Calculate coeff
f = DistBuilderExpCos2Dsing(Nc,gridObj.phi,Coeff_best);        % Build equil distribution
% plot(gridObj.phi,f)

% Initialize rho
rho = systemObj.numPart / ( systemObj.Lx .* systemObj.Ly) .* ...
    ones(systemObj.Nx,systemObj.Ny,systemObj.Nm);

% Map distribution to a homogeneous system
for i = 1:systemObj.Nm
    rho(:,:,i) = rho(:,:,i) .* f(i);
end

% systemObj.numPart / (systemObj.Lx .* systemObj.Lx);
% b = particleObj.lMaj^2 / pi;
% c =  systemObj.bc / b;
% f_reshape =  reshape(rho(17,17,:) / c  , 1, 32 );
% trapz_periodic(gridObj.phi,f)
% trapz_periodic(gridObj.phi,f_reshape)
% keyboard
% Normalize it
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(gridObj.y,trapz_periodic(gridObj.x,trapz_periodic(gridObj.phi,rho,3),2),1);
rho_eq = rho .* systemObj.numPart ./ CurrentNorm;
% Perturb it

% keyboard
end %end function