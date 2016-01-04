% Function: dRhoDriveCalcFT_ID
% 
% Description: Calculates the driving term contribution to drho/dt. This
% term should really be in the propagator. This is approximation that can
% be used by some AB stepping method.


function [dRhoDrive_FT] = dRhoDriveCalcFtId(rho,v,phi3D,kx3D,ky3D)

vx = v .* cos(phi3D);
vy = v .* sin(phi3D);
dRhoDrive_FT = - sqrt(-1) .* kx3D .* fftshift(fftn( rho .* vx ) ) + ... ;
    - sqrt(-1) .* ky3D .* fftshift(fftn( rho .* vy ) );

end

