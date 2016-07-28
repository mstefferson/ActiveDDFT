% PwDenPerturber2Drot.m
% Perturbs input rhp
% Perturbations of the form
%  \sum epsilon_k exp( i (k_x x + k_y y + k_phi phi)
% epsilon_k is an input parameter
function [rho] = PwDenPerturber2Drot(rho,ParamObj,GridObj,RhoInit)

[GridObj.y3D, GridObj.x3D, GridObj.phi3D] = ...
  meshgrid( GridObj.y, GridObj.x, GridObj.phi);
%% THIS WILL BREAK WITHOUT x3d. Need to simultaneously switch int rho method
% keyboard
MaxPerturb = RhoInit.WeightPos * ...
  (2*RhoInit.NumModesX)* (2*RhoInit.NumModesY)* (2*RhoInit.NumModesM);
if min(min(min(rho))) < MaxPerturb
  Coeff = RhoInit.WeightPos * min(min(min(rho))) / MaxPerturb;
else
  Coeff = RhoInit.WeightPos;
end

try
  % keyboard
  for i = -RhoInit.NumModesX:RhoInit.NumModesX
    for j = -RhoInit.NumModesY:RhoInit.NumModesY
      for k = -RhoInit.NumModesM: RhoInit.NumModesM
        %Change in x
        %                 keyboard
        if RhoInit.RandomAmp
          Coeff = (-1 + 2 * rand() ) * RhoInit.WeightPos;
        end
        rho = rho + ...
          ( Coeff + sqrt(-1) * Coeff ) .* (  ...
          exp(sqrt(-1) .* (...
          GridObj.kx(ParamObj.Nx/2+1+i) .* GridObj.x3D + ...
          GridObj.ky(ParamObj.Ny/2+1+j) .* GridObj.y3D + ...
          GridObj.km(ParamObj.Nm/2+1+k) .* GridObj.phi3D ) ) ) ;
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
CurrentNorm = trapz_periodic(GridObj.x,trapz_periodic(GridObj.y,trapz_periodic(GridObj.phi,rho,3),2),1);
rho = rho .* ParamObj.Norm ./ CurrentNorm;

end
