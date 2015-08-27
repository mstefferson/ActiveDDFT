% IntDen2DrotCalcPerturbEqSepPw.m
%
% Description: creates the initial density in 2 spatial directions and one
% rotation DoF. The initial density of the form
% rho = rho_equilibrium + \sum A_k exp( i k_x x ) + \sum A_k exp( i k_y y ) ...
%       + \sum A_k exp( i k_m phi )
% A_k is an input parameter

function [rho] = IntDenCalcPerturbEqSepPw2Drot(GridObj,ParamObj)

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

% Distribution stuff
Nc    = 10;            % Number of Coefficients

[Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,GridObj.phi,ParamObj.bc); % Calculate coeff
f = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        % Build equil distribution
% plot(GridObj.phi,f)

% Initialize rho
rho = ParamObj.Norm / (2 .* pi .* ParamObj.Lx .* ParamObj.Lx) .* ...
    ones(ParamObj.Nx,ParamObj.Ny,ParamObj.Nm);

% Map distribution to a homogeneous system
for i = 1:ParamObj.Nm
    rho(:,:,i) = rho(:,:,i) .* f(i);
end


% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(GridObj.y,trapz_periodic(GridObj.x,trapz_periodic(GridObj.phi,rho,3),2),1);
rho = rho .* ParamObj.Norm ./ CurrentNorm;
% Perturb it
[rho] = SepPwDenPerturber2Drot(rho,ParamObj,GridObj);

end %end function