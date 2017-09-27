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
n1 = systemObj.n1;
n2 = systemObj.n2;
n3 = systemObj.n3;

% Amplitudes
aX   = rhoInit.aXf .* systemObj.c;
aY   = rhoInit.aYf .* systemObj.c;
aPhi = rhoInit.aPhif .* systemObj.c;

% Build angular gaussian perturbation
centerPhi = gridObj.x3( rho(1,1,:) == max(rho(1,1,:) ) );
if length(centerPhi) > 1; centerPhi = centerPhi(1); end;
pPhi = gaussShift( gridObj.x3, systemObj.l3, aPhi, rhoInit.varPhi, centerPhi ); 
pPhi3 = repmat( reshape( pPhi, [1, 1, n3] ), [n1, n2, 1] );

% Build spatial gaussian perturbation (x)
if aX == 0
  pX3 = 1;
else
  pX = gaussShift( gridObj.x1, systemObj.l1, aX, rhoInit.varX , rhoInit.centerX); 
  pX3 = repmat( pX', [1, n2, n3] );
end

% Build spatial gaussian perturbation (y)
if aY == 0
  pY3 = 1;
else
  pY = gaussShift( gridObj.x2, systemObj.l2, aY, rhoInit.varY , rhoInit.centerY); 
  pY3 = repmat( pY, [n1, 1, n3] );
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
  normSqr = trapz_periodic( gridObj.x3, rho, 3 );
  normTemp = repmat( normSqr, [1, 1, n3 ] );
  rho = systemObj.numPart ./ normTemp .* rho;
else % Inhomogenous
  normTemp = trapz_periodic( gridObj.x1 , ...
    trapz_periodic( gridObj.x2, trapz_periodic( gridObj.x3, rho, 3 ), 2 ), 1 );
  rho = systemObj.numPart ./ normTemp .* rho;
end

end
