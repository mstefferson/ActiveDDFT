% PwDenPerturber2Drot.m
% Perturbs input rhp
% Perturbations of the form
%  \sum epsilon_k exp( i (k_x x + k_y y + k_phi phi)
% epsilon_k is an input parameter
function [rho] = PwDenPerturber2Drot(rho,ParamObj,GridObj)

% keyboard
MaxPerturb = ParamObj.WeightPos * ...
    (2*ParamObj.NumModesX)* (2*ParamObj.NumModesY)* (2*ParamObj.NumModesM);
if min(min(min(rho))) < MaxPerturb
Coeff = ParamObj.WeightPos * min(min(min(rho))) / MaxPerturb;
else
Coeff = ParamObj.WeightPos;
end

try
% keyboard
    for i = -ParamObj.NumModesX:ParamObj.NumModesX
        for j = -ParamObj.NumModesY:ParamObj.NumModesY
            for k = -ParamObj.NumModesM: ParamObj.NumModesM
                %Change in x
%                 keyboard
if ParamObj.Random
   Coeff = (-1 + 2 * rand() ) * ParamObj.WeightPos;
end
                rho = rho + ...
                    ( Coeff + sqrt(-1) * Coeff ) .* (  ...
                    exp(sqrt(-1) .* (...
                    GridObj.kx(ParamObj.Nx/2+1+i) .* GridObj.x3D + ...
                    GridObj.ky(ParamObj.Ny/2+1+j) .* GridObj.y3D + ...
                    GridObj.km(ParamObj.Nm/2+1+k) .* GridObj.phi3D ) ) ) ;
                
%             keyboard
%                 Perturb_FT = fftshift( fftn( ( Coeff + sqrt(-1) * Coeff ) .* (  ...
%                     exp(sqrt(-1) .* (...
%                     GridObj.kx(ParamObj.Nx/2+1+i) .* GridObj.x3D + ...
%                     GridObj.ky(ParamObj.Ny/2+1+j) .* GridObj.y3D + ...
%                     GridObj.km(ParamObj.Nm/2+1+k) .* GridObj.phi3D ) ) )  ) );
%                 Perturb_FT(ParamObj.Nx/2+1+i,ParamObj.Ny/2+1+j,ParamObj.Nm/2+1+k)
%                 keyboard
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
rho = rho .* ParamObj.Norm ./ CurrentNorm;

end