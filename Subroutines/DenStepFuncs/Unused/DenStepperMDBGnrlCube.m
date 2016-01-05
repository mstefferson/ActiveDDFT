function [rho_FT_next] = DenStepperMDBGnrlCube(Prop,rho_FT, GammaEx_FT,GammaEx_FT_prev,dt)

% Step using the hybrid AB. Everything is still in cube form.
rho_FT_next = Prop .* ( rho_FT ...
    + dt/2 * ( 3 * GammaEx_FT - Prop .* GammaEx_FT_prev ) );

end