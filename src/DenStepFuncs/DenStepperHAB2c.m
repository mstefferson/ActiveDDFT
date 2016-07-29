% Hybrid AB2

function [rho_FT_next] = DenStepperHAB2c(Prop,rho_FT, GammaEx_FT,GammaEx_FT_prev,dt)

% Step using the hybrid AB. Everything is still in cube form.
rho_FT_next = Prop .* ( rho_FT ...
    + dt * ( 3/2 * GammaEx_FT - 1/2 * Prop .* GammaEx_FT_prev ) );

end