% Function: dRhoDriveCalcFT_ID
% 
% Description: Calculates the driving term contribution to drho/dt. This
% term should really be in the propagator. This is approximation that can
% be used by some AB stepping method.


function [dRhoDrive_FT] = dRhoDriveCalcFtId(rho,v,...
  cosPhi3,sinPhi3,ikx3D,iky3D)

vx = v .* cosPhi3;
vy = v .* sinPhi3;
dRhoDrive_FT = - ikx3D .* fftshift(fftn( rho .* vx ) ) + ... ;
    - iky3D .* fftshift(fftn( rho .* vy ) );

end

