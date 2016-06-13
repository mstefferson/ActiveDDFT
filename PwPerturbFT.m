% PwPerturbFT.m
% Perturbs input rhp
% Perturbations of the form
%  \sum epsilon_k exp( i (k_x x + k_y y + k_phi phi)
% epsilon_k is an input parameter
% Does everything in Fourier space

function [rho] = PwPerturbFT(rho,ParamObj,RhoInit)

% N3 can be commonly used. Declare it
N3 = ParamObj.Nx * ParamObj.Ny * ParamObj.Nm;

% Perturb coeff is the weight times equilbrium
% concentration. Make sure it isn't too large
% Find isotropic density
IsoDen = ParamObj.Norm / (2 .* pi .* ParamObj.Lx .* ParamObj.Lx);

MaxPerturb = IsoDen * RhoInit.WeightPert * ...
  (2*RhoInit.NumModesX) * (2*RhoInit.NumModesY) * (2*RhoInit.NumModesM);

if min(min(min(rho))) < MaxPerturb
  CoeffMax = IsoDen .* ...
    RhoInit.WeightPert * min(min(min(rho))) / MaxPerturb;
else
  CoeffMax = IsoDen .* RhoInit.WeightPert; 
end

% If it's not random, set Coeff ourside the loop
% scale by N3 b/c of FT factor
if RhoInit.RandomAmp == 0
  Coeff = CoeffMax * N3 * ( 1 + sqrt(-1) );
else
  CoeffTemp = CoeffMax * N3;
end

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
          if RhoInit.RandomAmp
            Coeff = CoeffTemp .* ... 
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

% Fix negative rho if that happened.
[rho] = FixNegDenFnc(rho);

end
