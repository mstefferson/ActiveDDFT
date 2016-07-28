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

N     = 32;
d_phi = 2*pi/N;
Lx    = 2*pi;
Ly    = Lx;
dx    = Lx/N;
x     = 0:dx:Lx-dx;
phi   = 0:d_phi:2*pi - d_phi;
Nc    = 20;            % Number of Coefficients
bc    = 1.6;
Norm  = 100;

[Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,phi,bc); % Calculate coeff
f = DistBuilderExpCos2Dsing(Nc,phi,Coeff_best);        % Build equil distribution
% plot(GridObj.phi,f)
f_circ = circshift(f',8)';

% a = (1:4)'
% a_circ = circshift(a,1)'
figure
plot(phi,f,phi,f_circ)
% Initialize rho
rho = Norm / (Lx*Ly).* ones(N,N,N);

% Map distribution to a homogeneous system

rho(i,j,k) = rho(i,j,k) .* f_circ(k);
for i = 1:N
    for j = 1:N
        for k = 1:N
            rho(i,j,k) = rho(i,j,k) .* f(k);
        end
    end
end

rho2 = ones(N,N,N);

for i = 1:N
    for j = 1:N
        if mod( (i+j), 2 ) == 0
            rho2(i,j,:) = f;
        else
            rho2(i,j,:) = f_circ;
        end
    end
end
rho2 = rho2 * Norm / (Lx*Ly);
figure
DenRecObj.Density
plot( phi,reshape(rho2(2,1,:),1,32),phi, reshape(rho2(2,2,:),1,32) )
plot( phi,reshape(rho(2,1,:),1,32),phi, reshape(rho(2,2,:),1,32) )
plot( phi,reshape(DenRecObj.Density_rec(2,1,:,1),1,32),...
    phi, reshape(DenRecObj.Density_rec(2,2,:,1),1,32) )
% ParamObj.Norm / (ParamObj.Lx .* ParamObj.Lx);
% b = ParamObj.L_rod^2 / pi;
% c =  ParamObj.bc / b;
% f_reshape =  reshape(rho(17,17,:) / c  , 1, 32 );
% trapz_periodic(GridObj.phi,f)
% trapz_periodic(GridObj.phi,f_reshape)
% keyboard
% Normalize it
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(GridObj.y,trapz_periodic(GridObj.x,trapz_periodic(GridObj.phi,rho,3),2),1);
rho_eq = rho .* ParamObj.Norm ./ CurrentNorm;
% Perturb it

% keyboard
