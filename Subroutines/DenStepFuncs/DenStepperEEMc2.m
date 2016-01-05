% Currently causes things to blow up
% Exponential Euler Method 2

function [rho_FT_next] = ...
    DenStepperEEMc2( Prop, GamProp, rho_FT, Gamma_FT, Gamma_FTprev)

%Take the first step in k-space using Euler. Everything is a cube.
%All this is being done element by element.

rho_FT_next = Prop .* rho_FT + GamProp .* ( 3 * Gamma_FT - Gamma_FTprev );

end