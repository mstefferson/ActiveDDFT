% PwDenPerturber2Drot.m
% Perturbs input rhp
% Perturbations of the form
%  \sum epsilon_k exp( i (k_x x + k_y y + k_phi phi)
% epsilon_k is an input parameter
function [rho] = PwDenPerturber2Drot(rho,ParamObj,gridObj,rhoInit)

[gridObj.x23D, gridObj.x13D, gridObj.x33D] = ...
  meshgrid( gridObj.x2, gridObj.x1, gridObj.x3);
%% THIS WILL BREAK WITHOUT x3d. Need to simultaneously switch int rho method
% keyboard
MaxPerturb = rhoInit.WeightPos * ...
  (2*rhoInit.NumModesX)* (2*rhoInit.NumModesY)* (2*rhoInit.NumModesM);
if min(min(min(rho))) < MaxPerturb
  Coeff = rhoInit.WeightPos * min(min(min(rho))) / MaxPerturb;
else
  Coeff = rhoInit.WeightPos;
end

try
  % keyboard
  for i = -rhoInit.NumModesX:rhoInit.NumModesX
    for j = -rhoInit.NumModesY:rhoInit.NumModesY
      for k = -rhoInit.NumModesM: rhoInit.NumModesM
        %Change in x
        %                 keyboard
        if rhoInit.RandomAmp
          Coeff = (-1 + 2 * rand() ) * rhoInit.WeightPos;
        end
        rho = rho + ...
          ( Coeff + sqrt(-1) * Coeff ) .* (  ...
          exp(sqrt(-1) .* (...
          gridObj.k1(systemObj.n1/2+1+i) .* gridObj.x13D + ...
          gridObj.k2(systemObj.n2/2+1+j) .* gridObj.x23D + ...
          gridObj.km(systemObj.n3/2+1+k) .* gridObj.x33D ) ) ) ;
      end
    end
  end % End loop over modes
  
catch err
  
  fprintf('%s', err.getReport('extended', 'hyperlinks','off')) ;
  fprintf('%s', err.getReport('extended')) ;
  keyboard
end
% Take real part
rho = real(rho);

% Fix negative rho if that happened.
[rho] = FixNegDenFnc(rho);

% Normalize it again just in case
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(gridObj.x1,trapz_periodic(gridObj.x2,trapz_periodic(gridObj.x3,rho,3),2),1);
rho = rho .* systemObj.numPart ./ CurrentNorm;

end
