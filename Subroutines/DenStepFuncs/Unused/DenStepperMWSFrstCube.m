function [rho_FT_next] = ...
    DenStepperMWSFrstCube( Prop, rho_FT, GammaEx_FT,dt)

%Take the first step in k-space using a hybrid Euler. Everything is a cube.
%All this is being done element by element.

rho_FT_next = Prop .* (rho_FT + dt/2 * GammaEx_FT) + dt/2 * GammaEx_FT;

end