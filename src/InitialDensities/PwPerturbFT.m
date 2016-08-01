% PwPerturbFT.m
% Perturbs input rhp
% Perturbations of the form
%  \sum epsilon_k exp( i (k_x x + k_y y + k_phi phi)
% epsilon_k is an input parameter
% Does everything in Fourier space

function [rho] = PwPerturbFT(rho,systemObj,rhoInit)
% N3 can be commonly used. Declare it
N3 = systemObj.Nx * systemObj.Ny * systemObj.Nm;
% Perturb coeff is the weight times equilbrium
% concentration. Make sure it isn't too large
% Find isotropic density
IsoDen = systemObj.numPart / (systemObj.Lphi .* systemObj.Lx .* systemObj.Lx);
MaxPerturb = IsoDen * rhoInit.WeightPert * ...
  (2*rhoInit.NumModesX) * (2*rhoInit.NumModesY) * (2*rhoInit.NumModesM);
if min(min(min(rho))) < MaxPerturb
  CoeffMax = IsoDen .* ...
    rhoInit.WeightPert * min(min(min(rho))) / MaxPerturb;
else
  CoeffMax = IsoDen .* rhoInit.WeightPert; 
end
% If it's not random, set Coeff ourside the loop
% scale by N3 b/c of FT factor
if rhoInit.RandomAmp == 0
  Coeff = CoeffMax * N3 * ( 1 + sqrt(-1) );
else
  CoeffTemp = CoeffMax * N3;
end
% Handle perturbations in Fourier space
rhoFT = fftshift( fftn( rho ) );
kx0   = systemObj.Nx / 2 + 1;
ky0   = systemObj.Ny / 2 + 1;
km0   = systemObj.Nm / 2 + 1;
try
  % Loop over perturbations
  for ii = 0:rhoInit.NumModesX
    for jj = 0:rhoInit.NumModesY
      for kk = 0: rhoInit.NumModesM
        if ii ~= 0 || jj ~=0 || kk ~= 0
          if rhoInit.RandomAmp
            Coeff = CoeffTemp .* ... 
              ( (-1 + 2 * rand() )  + (-1 + 2 * rand() ) * sqrt(-1) ); 
          end
          rhoFT(kx0 + ii, ky0 + jj, km0 + kk) =  ...
            rhoFT(kx0 + ii, ky0 + jj, km0 + kk) + Coeff;
          rhoFT(kx0 - ii, ky0 - jj, km0 - kk) =  ...
            rhoFT(kx0 - ii, ky0 - jj, km0 - kk) + conj( Coeff );
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

% Fix negative rho if that happened.
[rho] = FixNegDenFnc(rho);

end
