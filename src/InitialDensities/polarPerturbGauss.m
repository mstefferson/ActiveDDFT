% [rho] = polarPerturbGauss( rho, systemObj, rhoInit, gridObj )
% Description: Add gaussian perturbation in phi and spatially depending on given
% amplitudes.
%
% Key:
% 4: Polar perturbation. Homogenous concentration
% 5: Polar perturbation. Inhomogenous concentration

% Spatial perurbations are added if they have a non-zero amplitude

function [rho] = polarPerturbGauss( rho, systemObj, rhoInit, gridObj )

% Commonly used variables
Nx = systemObj.Nx;
Ny = systemObj.Ny;
Nm = systemObj.Nm;

% Amplitudes
aX   = rhoInit.aXf .* systemObj.c;
aY   = rhoInit.aYf .* systemObj.c;
aPhi = rhoInit.aPhif .* systemObj.c;

% Build angular gaussian perturbation
centerPhi = gridObj.phi( rho(1,1,:) == max(rho(1,1,:) ) );
if length(centerPhi) > 1; centerPhi = centerPhi(1); end;
pPhi = gaussShift( gridObj.phi, systemObj.Lphi, aPhi, rhoInit.varPhi, centerPhi ); 
pPhi3 = repmat( reshape( pPhi, [1, 1, Nm] ), [Nx, Ny, 1] );

% Build spatial gaussian perturbation (x)
if aX == 0
  pX3 = 1;
else
  pX = gaussShift( gridObj.x, systemObj.Lx, aX, rhoInit.varX , rhoInit.centerX); 
  pX3 = repmat( pX', [1, Ny, Nm] );
end

% Build spatial gaussian perturbation (y)
if aY == 0
  pY3 = 1;
else
  pY = gaussShift( gridObj.y, systemObj.Ly, aY, rhoInit.varY , rhoInit.centerY); 
  pY3 = repmat( pY, [Nx, 1, Nm] );
end

% Total perturbation
pTot = pX3 .* pY3 .* pPhi3;
rho = rho + pTot;
% Fix negative rho
minRho = min( min( min( rho ) ) );
if minRho < 0;
  rho = rho - minRho;
end

% Normalize base on homogenous or not
if rhoInit.IntCond == 4; % Homo
  normSqr = trapz_periodic( gridObj.phi, rho, 3 );
  normTemp = repmat( normSqr, [1, 1, Nm ] );
  rho = systemObj.numPart ./ normTemp .* rho;
else % Inhomogenous
  normTemp = trapz_periodic( gridObj.x , ...
    trapz_periodic( gridObj.y, trapz_periodic( gridObj.phi, rho, 3 ), 2 ), 1 );
  rho = systemObj.numPart ./ normTemp .* rho;
end

end
