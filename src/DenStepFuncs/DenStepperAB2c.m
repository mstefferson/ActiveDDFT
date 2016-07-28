
% Adams Bashforth 2
% Ignores integrating factor.

function [rho_FT_next] = ...
    DenStepperAB2c( Prop, rho_FT, GammaEx_FT,GammaEx_FTprev, dt)

%Take the first step in k-space using Euler. Everything is a cube.
%All this is being done element by element.

rho_FT_next = Prop .* rho_FT + dt/2 * ( 3 * GammaEx_FT -  GammaEx_FTprev );

end
