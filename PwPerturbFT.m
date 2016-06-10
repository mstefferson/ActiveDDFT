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
  CoeffMax = RhoInit.WeightPos * min(min(min(rho))) / MaxPerturb;
else
  CoeffMax = RhoInit.WeightPos;
end

% If it's not random, set Coeff ourside the loop
Coeff = CoeffMax * ( 1 + sqrt(-1) );
% Handle perturbations in Fourier space
rhoFT = fftshift( fftn( rho ) );
kx0   = ParamObj.Nx / 2 + 1;
ky0   = ParamObj.Ny / 2 + 1;
km0   = ParamObj.Nm / 2 + 1;

try

  
  % Loop over perturbations
  for ii = 0:RhoInit.NumModesX
    for jj = 0:RhoInit.NumModesY
      for kk = 0: RhoInit.NumModesM
        
        if ii ~= 0 || jj ~=0 || kk ~= 0
          fprintf('Perturbing (%d,%d,%d)',ii,jj,kk);
          %rhoFT(ii,jj,kk) =  Coeff;
          if RhoInit.RandomAmp
          Coeff = ( CoeffMax ) / N3 * ...
            ( (-1 + 2 * rand() )  + (-1 + 2 * rand() ) * sqrt(-1) ); 
          end
          rhoFT(kx0 + ii, ky0 + jj, km0 + kk) =  Coeff;
          rhoFT(kx0 - ii, ky0 - jj, km0 - kk) =  conj( Coeff );

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
end
