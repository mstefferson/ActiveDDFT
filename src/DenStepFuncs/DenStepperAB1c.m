% Adams Bashforth 1 (Euler)
% Ignores integrating factor.

function [rho_FT_next] = ...
    DenStepperAB1c( Prop, rho_FT, GammaEx_FT,dt)

%Take the first step in k-space using Euler. Everything is a cube.
%All this is being done element by element.

rho_FT_next = Prop .* rho_FT + dt * GammaEx_FT;

end