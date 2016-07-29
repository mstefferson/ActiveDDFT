% PwDenPerturber2Drot.m
% Perturbs input rhp
% Perturbations of the form
%  \sum epsilon_k exp( i (k_x x + k_y y + k_phi phi)
% epsilon_k is an input parameter
function [rho] = PwDenPerturber2Drot(rho,ParamObj,gridObj,rhoInit)

[gridObj.y3D, gridObj.x3D, gridObj.phi3D] = ...
  meshgrid( gridObj.y, gridObj.x, gridObj.phi);
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
          gridObj.kx(systemObj.Nx/2+1+i) .* gridObj.x3D + ...
          gridObj.ky(systemObj.Ny/2+1+j) .* gridObj.y3D + ...
          gridObj.km(systemObj.Nm/2+1+k) .* gridObj.phi3D ) ) ) ;
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
CurrentNorm = trapz_periodic(gridObj.x,trapz_periodic(gridObj.y,trapz_periodic(gridObj.phi,rho,3),2),1);
rho = rho .* systemObj.numPart ./ CurrentNorm;

end
