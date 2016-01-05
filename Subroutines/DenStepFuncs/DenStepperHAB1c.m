% Hybrid AB 1

function [rho_FT_next] = ...
    DenStepperHAB1c( Prop, rho_FT, GammaEx_FT,dt)

%Take the first step in k-space using Euler. Everything is a cube.
%All this is being done element by element.

rho_FT_next = Prop .* (rho_FT + dt * GammaEx_FT);

end