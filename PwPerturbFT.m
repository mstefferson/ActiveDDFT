% PwDenPerturber2Drot.m
% Perturbs input rhp
% Perturbations of the form
%  \sum epsilon_k exp( i (k_x x + k_y y + k_phi phi)
% epsilon_k is an input parameter
function [rho] = PwPerturbFT(rho,ParamObj,GridObj,RhoInit)

N3 = ParamObj.Nx * ParamObj.Ny * ParamObj.Nm;
MaxPerturb = RhoInit.WeightPos * ...
  (2*RhoInit.NumModesX)* (2*RhoInit.NumModesY)* (2*RhoInit.NumModesM);
if min(min(min(rho))) < MaxPerturb
  Coeff = RhoInit.WeightPos * min(min(min(rho))) / MaxPerturb;
else
  Coeff = RhoInit.WeightPos;
end

% Handle perturbations in Fourier space
rhoFT = fftshift( fftn( rho ) );
kx0   = ParamObj.Nx / 2 + 1;
ky0   = ParamObj.Ny / 2 + 1;
km0   = ParamObj.Nm / 2 + 1;

try
  % keyboard
  % Loop over perturbations
  for ii = -RhoInit.NumModesX:RhoInit.NumModesX
    for jj = -RhoInit.NumModesY:RhoInit.NumModesY
      for kk = -RhoInit.NumModesM: RhoInit.NumModesM
        
        if ii ~= 0 || jj ~=0 || kk ~= 0
          %rhoFT(ii,jj,kk) =  Coeff;
          rhoFT(kx0 + ii, ky0 + jj, km0 + kk) =  Coeff / N3;
        end
        
      end
    end
  end % End loop over modes
  
catch err
  
  fprintf('%s', err.getReport('extended', 'hyperlinks','off')) ;
  fprintf('%s', err.getReport('ex:tended')) ;
  keyboard
end
% Inverse transform and Take real part
rho = real( ifftn( ifftshift( rhoFT ) ) );

%     figure
%     surf( rho(:,:,17) )
%
% keyboard
% Fix negative rho if that happened.
[rho] = FixNegDenFnc(rho);

% Normalize it again just in case
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(GridObj.x,trapz_periodic(GridObj.y,trapz_periodic(GridObj.phi,rho,3),2),1);
disp( CurrentNorm ); disp( ParamObj.Norm );
end
