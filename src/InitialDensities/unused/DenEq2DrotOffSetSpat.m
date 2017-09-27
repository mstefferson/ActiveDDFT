% IntDen2DrotCalcPerturbEqPw.m 
% Description: creates the initial density in 2 spatial directions and one
% rotation DoF. The initial density of the form
% rho = rho_equilbrium + \sum A_k exp( i (k_x x + k_y y + k_phi phi)
% A_k is an input parameter


function [rho] = DenEq2DrotOffSetSpat(gridObj,ParamObj)

% Add in some slight deviation from the equilbrium density at specific modes.
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

[Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,gridObj.x3,systemObj.bc); % Calculate coeff
f = DistBuilderExpCos2Dsing(Nc,gridObj.x3,Coeff_best);        % Build equil distribution
% Get a new dist with peaks shifted by pi/2
f_circ = circshift(f',systemObj.n3/4)';

% Initialize rho
rho =    ones(systemObj.n1,systemObj.n2,systemObj.n3);

% Map distribution to a homogeneous system
for i = 1:systemObj.n1
    for j = 1:systemObj.n2
        if mod( (i+j), 2 ) == 0
            rho(i,j,:) = f;
        else
            rho(i,j,:) = f_circ;
        end
    end
end

% Normalize it
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(gridObj.x2,trapz_periodic(gridObj.x1,trapz_periodic(gridObj.x3,rho,3),2),1);
rho = rho .* systemObj.numPart ./ CurrentNorm;
% Perturb it

% keyboard
end %end function