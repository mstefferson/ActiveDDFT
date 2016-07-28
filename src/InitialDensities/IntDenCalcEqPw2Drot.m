% IntDen2DrotCalcPerturbEqPw.m 
% Description: creates the initial density in 2 spatial directions and one
% rotation DoF. The initial density of the form
% rho = rho_equilbrium + \sum A_k exp( i (k_x x + k_y y + k_phi phi)
% A_k is an input parameter


function [rho] = IntDenCalcEqPw2Drot(ParamObj, RhoInit, x, y, phi)

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
rho = ParamObj.c .* ...
    ones(ParamObj.Nx,ParamObj.Ny,ParamObj.Nm);

% Map distribution to a homogeneous system
for i = 1:ParamObj.Nm
    rho(:,:,i) = rho(:,:,i) .* RhoInit.feq(i);
end

% ParamObj.Norm / (ParamObj.Lx .* ParamObj.Lx);
% b = ParamObj.L_rod^2 / pi;
% c =  ParamObj.bc / b;
% f_reshape =  reshape(rho(17,17,:) / c  , 1, 32 );
% trapz_periodic(phi,f)
% trapz_periodic(phi,f_reshape)
% keyboard
% Normalize it
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(x,trapz_periodic(y,trapz_periodic(phi,rho,3),2),1);
rho = rho .* ParamObj.Norm ./ CurrentNorm;

% Perturb it
%[rho] = PwDenPerturber2Drot(rho_eq,ParamObj,RhoInit);
[rho] = PwPerturbFT(rho,ParamObj,RhoInit);

end %end function
