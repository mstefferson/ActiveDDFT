% PwDenPerturber2Drot.m
% Perturbs input rhp
% Perturbations of the form
% \sum epsilon_k exp( i k_x x ) + \sum epsilon_k exp( i k_y y ) ...
%       + \sum epsilon_k exp( i k_m phi )

function [rho] = SepPwDenPerturber2Drot(rho,ParamObj,gridObj,rhoInit)

% NEED TO EDIT. This will break with changes to gridObj
%Perturbation coeff
Coeff = rhoInit.WeightPos;
%Change in x
for i = -rhoInit.NumModesX:rhoInit.NumModesX
  
  
  rho = rho + (Coeff + sqrt(-1) .* Coeff) ...
    .* exp( sqrt(-1) .* gridObj.k1(systemObj.n1/2+1+i) .* gridObj.x13D ) ; ...
end

%Change in y
for i = -rhoInit.NumModesY:rhoInit.NumModesY
  
  rho = rho + (Coeff + sqrt(-1) .* Coeff) ...
    .* exp( sqrt(-1) .* gridObj.k2(systemObj.n2/2+1+i) .* gridObj.x23D ) ; ...
end

Coeff = rhoInit.WeightAng;
%Change in phi
for i = -rhoInit.NumModesM:rhoInit.NumModesM
  
  rho = rho + (Coeff + sqrt(-1) .* Coeff) ...
    .* exp( sqrt(-1) .* gridObj.km(systemObj.n3/2+1+i) .* gridObj.x33D ) ; ...
end

% Take real part
rho = real(rho);

% Fix negative rho if that happened.
[rho] = FixNegDenFnc(rho);

% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(gridObj.x2,trapz_periodic(gridObj.x1,trapz_periodic(gridObj.x3,rho,3),2),1);
rho = real(rho .* systemObj.numPart ./ CurrentNorm);

end
